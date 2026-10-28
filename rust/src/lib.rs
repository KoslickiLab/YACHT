//! Native core of YACHT: a compact, memory-mapped index of a trained reference database and the
//! queries `yacht run` needs (which genomes overlap a sample, and how many of each genome's exclusive
//! k-mers are present in the sample).

pub mod builder;
pub mod format;
pub mod query;
pub mod sig;
pub mod tsv;
pub mod xlsx;

#[cfg(feature = "python")]
mod python {
    use crate::builder;
    use crate::tsv;
    use crate::xlsx;
    use crate::query::{Index, QueryResult};
    use numpy::{PyArray1, PyReadonlyArray1};
    use pyo3::exceptions::{PyRuntimeError, PyValueError};
    use pyo3::prelude::*;
    use pyo3::types::{PyBytes, PyDict};
    use std::path::{Path, PathBuf};

    /// Memory-mapped YACHT run index.
    #[pyclass(name = "YachtIndex", frozen, module = "yacht._yacht_core")]
    struct PyIndex {
        inner: Index,
    }

    /// Result of looking up a sample's hashes in a YachtIndex.
    #[pyclass(name = "IndexQuery", frozen, module = "yacht._yacht_core")]
    struct PyQuery {
        inner: QueryResult,
    }

    #[pymethods]
    impl PyQuery {
        /// Number of sample hashes present in the index.
        #[getter]
        fn n_found(&self) -> u64 {
            self.inner.n_found
        }
    }

    #[pymethods]
    impl PyIndex {
        #[new]
        fn new(path: PathBuf) -> PyResult<Self> {
            Index::open(&path).map(|inner| PyIndex { inner }).map_err(PyValueError::new_err)
        }

        #[getter]
        fn format_version(&self) -> u32 {
            self.inner.hdr.version
        }

        #[getter]
        fn ksize(&self) -> u64 {
            self.inner.hdr.ksize
        }

        #[getter]
        fn scaled(&self) -> u64 {
            self.inner.hdr.scaled
        }

        #[getter]
        fn n_rows(&self) -> u64 {
            self.inner.hdr.n_rows
        }

        #[getter]
        fn n_keys(&self) -> u64 {
            self.inner.hdr.n_keys
        }

        #[getter]
        fn n_sets(&self) -> u64 {
            self.inner.hdr.n_sets
        }

        /// SHA-256 digest of the manifest this index was built for.
        #[getter]
        fn digest<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
            PyBytes::new(py, &self.inner.hdr.digest)
        }

        /// Number of hashes of each genome sketch, in manifest row order.
        fn sketch_sizes<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
            PyArray1::from_slice(py, self.inner.sizes())
        }

        /// Signature names (as reported by multisearch) of the given rows.
        fn names(&self, rows: PyReadonlyArray1<'_, u32>) -> PyResult<Vec<String>> {
            let n = self.inner.hdr.n_rows as usize;
            rows.as_slice()?
                .iter()
                .map(|&r| {
                    if (r as usize) < n {
                        self.inner.name(r as usize).map_err(PyRuntimeError::new_err)
                    } else {
                        Err(PyValueError::new_err("row out of range"))
                    }
                })
                .collect()
        }

        /// Phase A: look up the sample hashes (uint64 array).
        fn query(&self, hashes: PyReadonlyArray1<'_, u64>) -> PyResult<PyQuery> {
            Ok(PyQuery {
                inner: self.inner.query(hashes.as_slice()?),
            })
        }

        /// Rows sharing at least one hash with the sample (sorted uint32 array).
        fn marked_rows<'py>(&self, py: Python<'py>, q: &PyQuery) -> Bound<'py, PyArray1<u32>> {
            PyArray1::from_vec(py, self.inner.marked_rows(&q.inner))
        }

        /// Phase B: (excl, matches) uint64 arrays of length n_rows for the nontrivial-row mask.
        fn exclusive_counts<'py>(
            &self,
            py: Python<'py>,
            q: &PyQuery,
            mask: PyReadonlyArray1<'_, bool>,
        ) -> PyResult<(Bound<'py, PyArray1<u64>>, Bound<'py, PyArray1<u64>>)> {
            let (excl, matches) = self
                .inner
                .exclusive_counts(&q.inner, mask.as_slice()?)
                .map_err(PyValueError::new_err)?;
            Ok((PyArray1::from_vec(py, excl), PyArray1::from_vec(py, matches)))
        }
    }

    /// Build a run index from genome signature files (one per manifest row, in manifest order).
    #[pyfunction]
    fn build_index<'py>(
        py: Python<'py>,
        sig_paths: Vec<PathBuf>,
        ksize: u64,
        scaled: u64,
        digest: &[u8],
        out_path: PathBuf,
    ) -> PyResult<Bound<'py, PyDict>> {
        let digest: [u8; 32] = digest
            .try_into()
            .map_err(|_| PyValueError::new_err("digest must be 32 bytes"))?;
        let stats = builder::build_index(&sig_paths, ksize, scaled, digest, Path::new(&out_path))
            .map_err(PyValueError::new_err)?;
        let d = PyDict::new(py);
        d.set_item("n_rows", stats.n_rows)?;
        d.set_item("n_hashes", stats.n_hashes)?;
        d.set_item("n_keys", stats.n_keys)?;
        d.set_item("n_unique_keys", stats.n_unique_keys)?;
        d.set_item("n_sets", stats.n_sets)?;
        d.set_item("n_members", stats.n_members)?;
        d.set_item("file_len", stats.file_len)?;
        Ok(d)
    }

    /// Column data extracted from Python, kept alive while the Rust columns borrow from it.
    enum Held<'py> {
        I(PyReadonlyArray1<'py, i64>),
        F(PyReadonlyArray1<'py, f64>),
        B(PyReadonlyArray1<'py, bool>),
        S(Vec<Option<String>>),
    }

    /// Extract (kind, values) column specs: kind "i" (int64 array), "f" (float64 array), "b" (bool array)
    /// or "s" (list of str or None).
    fn hold_columns<'py>(columns: &[(String, Bound<'py, PyAny>)]) -> PyResult<Vec<Held<'py>>> {
        let mut held: Vec<Held> = Vec::with_capacity(columns.len());
        for (kind, values) in columns {
            held.push(match kind.as_str() {
                "i" => Held::I(values.extract()?),
                "f" => Held::F(values.extract()?),
                "b" => Held::B(values.extract()?),
                "s" => Held::S(values.extract()?),
                _ => return Err(PyValueError::new_err(format!("unknown column kind {}", kind))),
            });
        }
        Ok(held)
    }

    fn as_columns<'a>(held: &'a mut [Held<'_>]) -> PyResult<Vec<tsv::Column<'a>>> {
        let mut cols: Vec<tsv::Column> = Vec::with_capacity(held.len());
        for h in held.iter_mut() {
            cols.push(match h {
                Held::I(a) => tsv::Column::Int(a.as_slice()?),
                Held::F(a) => tsv::Column::Float(a.as_slice()?),
                Held::B(a) => tsv::Column::Bool(a.as_slice()?),
                Held::S(v) => tsv::Column::Text(std::mem::take(v)),
            });
        }
        Ok(cols)
    }

    /// Write a table as tab-separated text exactly like `pandas.DataFrame.to_csv(path, sep="\t", index=False)`.
    /// `columns` is a list of (kind, values), see hold_columns.
    #[pyfunction]
    fn write_tsv(
        path: PathBuf,
        header: Vec<String>,
        columns: Vec<(String, Bound<'_, PyAny>)>,
        line_terminator: &str,
    ) -> PyResult<()> {
        let mut held = hold_columns(&columns)?;
        for h in &held {
            if let Held::S(v) = h {
                if v.iter().flatten().any(|s| !tsv::text_supported(s)) {
                    return Err(PyValueError::new_err("text with line breaks is not supported"));
                }
            }
        }
        if header.iter().any(|s| !tsv::text_supported(s)) {
            return Err(PyValueError::new_err("text with line breaks is not supported"));
        }
        let cols = as_columns(&mut held)?;
        tsv::write_tsv(&path, &header, &cols, line_terminator).map_err(PyRuntimeError::new_err)
    }

    /// Write an .xlsx workbook from the parts of a template workbook. `parts` is a list of
    /// (name, head, tail, sheet) where the data rows of sheet number `sheet` (-1: none) are written between
    /// head and tail; `sheets` holds the columns of each sheet, see hold_columns.
    #[pyfunction]
    fn write_xlsx(
        path: PathBuf,
        parts: Vec<(String, Vec<u8>, Vec<u8>, i64)>,
        sheets: Vec<Vec<(String, Bound<'_, PyAny>)>>,
    ) -> PyResult<()> {
        let mut held: Vec<Vec<Held>> = sheets.iter().map(|c| hold_columns(c)).collect::<PyResult<_>>()?;
        let mut cols: Vec<Vec<tsv::Column>> = Vec::with_capacity(held.len());
        for h in held.iter_mut() {
            cols.push(as_columns(h)?);
        }
        let parts: Vec<xlsx::Part> = parts
            .into_iter()
            .map(|(name, head, tail, sheet)| xlsx::Part {
                name,
                head,
                tail,
                sheet: usize::try_from(sheet).ok(),
            })
            .collect();
        xlsx::write_xlsx(&path, &parts, &cols).map_err(PyRuntimeError::new_err)
    }

    /// Load the sketch of a sample signature file (the bytes of one `.sig`/`.sig.gz` entry of a sample zip):
    /// the DNA sketch with k-mer size `ksize` whose sourmash md5sum is `md5`. Returns
    /// (hashes uint64 array, abundances uint64 array or None, scaled, name), or None if the file is not a plain,
    /// unambiguous scaled sketch (the caller then loads it with sourmash).
    #[pyfunction]
    fn load_sample_sketch<'py>(
        py: Python<'py>,
        data: &[u8],
        ksize: u64,
        md5: &str,
    ) -> PyResult<Option<(Bound<'py, PyArray1<u64>>, Option<Bound<'py, PyArray1<u64>>>, u64, String)>> {
        let parsed = crate::sig::parse_sample_sketch(data, ksize, md5).map_err(PyValueError::new_err)?;
        Ok(parsed.map(|s| {
            (
                PyArray1::from_vec(py, s.hashes),
                s.abundances.map(|a| PyArray1::from_vec(py, a)),
                s.scaled,
                s.name,
            )
        }))
    }

    /// Python "%.16g" formatting of each value of a float64 array, as written by write_xlsx. For testing.
    #[pyfunction]
    fn format_floats_g16(values: PyReadonlyArray1<'_, f64>) -> PyResult<Vec<String>> {
        Ok(values
            .as_slice()?
            .iter()
            .map(|&x| {
                let mut v = Vec::new();
                xlsx::push_g16(&mut v, x);
                String::from_utf8(v).unwrap()
            })
            .collect())
    }

    /// Python repr() of each value of a float64 array (NaN -> ""), as written by write_tsv. For testing.
    #[pyfunction]
    fn format_floats(values: PyReadonlyArray1<'_, f64>) -> PyResult<Vec<String>> {
        Ok(values
            .as_slice()?
            .iter()
            .map(|&x| {
                let mut v = Vec::new();
                if !x.is_nan() {
                    tsv::push_float_repr(&mut v, x);
                }
                String::from_utf8(v).unwrap()
            })
            .collect())
    }

    #[pymodule]
    fn _yacht_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<PyIndex>()?;
        m.add_class::<PyQuery>()?;
        m.add_function(wrap_pyfunction!(build_index, m)?)?;
        m.add_function(wrap_pyfunction!(write_tsv, m)?)?;
        m.add_function(wrap_pyfunction!(format_floats, m)?)?;
        m.add_function(wrap_pyfunction!(write_xlsx, m)?)?;
        m.add_function(wrap_pyfunction!(format_floats_g16, m)?)?;
        m.add_function(wrap_pyfunction!(load_sample_sketch, m)?)?;
        m.add("FORMAT_VERSION", crate::format::FORMAT_VERSION)?;
        Ok(())
    }
}
