//! Fast writer for `result.xlsx`.
//!
//! The Python caller lets pandas/openpyxl write a template workbook with the same sheets but only the header
//! rows. This module writes the final .xlsx: every part of the template is copied as is, except that the data
//! rows are inserted into each worksheet. Cells are encoded exactly as openpyxl encodes the values that
//! `DataFrame.to_excel` gives it, so reading the file back yields the same cell values and types.

use crate::tsv::Column;
use flate2::write::DeflateEncoder;
use flate2::{Compression, Crc};
use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

/// One part (file) of the workbook. For worksheets, the data rows of `sheet` are written between
/// `head` and `tail`.
pub struct Part {
    pub name: String,
    pub head: Vec<u8>,
    pub tail: Vec<u8>,
    pub sheet: Option<usize>,
}

/// Append Python's `"%.16g" % x` (what openpyxl writes for numbers) for a finite double.
pub fn push_g16(out: &mut Vec<u8>, x: f64) {
    // 16 significant digits, correctly rounded (ties to even, like Python's formatting)
    let s = format!("{:.15e}", x);
    let (mant, exp) = s.split_once('e').expect("float formatted without exponent");
    let exp: i32 = exp.parse().expect("bad exponent");
    let (neg, mant) = match mant.strip_prefix('-') {
        Some(m) => (true, m),
        None => (false, mant),
    };
    let mut digits: Vec<u8> = mant.bytes().filter(|&c| c != b'.').collect();
    while digits.len() > 1 && digits.last() == Some(&b'0') {
        digits.pop();
    }
    if neg {
        out.push(b'-');
    }
    if digits == [b'0'] {
        out.push(b'0');
        return;
    }
    if exp < -4 || exp >= 16 {
        out.push(digits[0]);
        if digits.len() > 1 {
            out.push(b'.');
            out.extend_from_slice(&digits[1..]);
        }
        out.push(b'e');
        out.push(if exp < 0 { b'-' } else { b'+' });
        let a = exp.unsigned_abs();
        if a < 10 {
            out.push(b'0');
        }
        out.extend_from_slice(a.to_string().as_bytes());
    } else if exp < 0 {
        out.extend_from_slice(b"0.");
        out.extend(std::iter::repeat(b'0').take((-exp - 1) as usize));
        out.extend_from_slice(&digits);
    } else {
        let int_len = exp as usize + 1;
        if digits.len() <= int_len {
            out.extend_from_slice(&digits);
            out.extend(std::iter::repeat(b'0').take(int_len - digits.len()));
        } else {
            out.extend_from_slice(&digits[..int_len]);
            out.push(b'.');
            out.extend_from_slice(&digits[int_len..]);
        }
    }
}

fn push_escaped(out: &mut Vec<u8>, s: &str) {
    for c in s.bytes() {
        match c {
            b'&' => out.extend_from_slice(b"&amp;"),
            b'<' => out.extend_from_slice(b"&lt;"),
            b'>' => out.extend_from_slice(b"&gt;"),
            _ => out.push(c),
        }
    }
}

fn col_name(mut c: usize) -> Vec<u8> {
    let mut v = Vec::new();
    c += 1;
    while c > 0 {
        let r = (c - 1) % 26;
        v.push(b'A' + r as u8);
        c = (c - 1) / 26;
    }
    v.reverse();
    v
}

fn push_inline_str(out: &mut Vec<u8>, cref: &[u8], s: &str) {
    out.extend_from_slice(b"<c r=\"");
    out.extend_from_slice(cref);
    if s.is_empty() {
        out.extend_from_slice(b"\" t=\"inlineStr\" />");
        return;
    }
    out.extend_from_slice(b"\" t=\"inlineStr\"><is><t");
    let stripped = s.trim();
    if !stripped.is_empty() && stripped.len() != s.len() {
        out.extend_from_slice(b" xml:space=\"preserve\"");
    }
    out.push(b'>');
    push_escaped(out, s);
    out.extend_from_slice(b"</t></is></c>");
}

fn push_number(out: &mut Vec<u8>, cref: &[u8], x: f64) {
    out.extend_from_slice(b"<c r=\"");
    out.extend_from_slice(cref);
    out.extend_from_slice(b"\" t=\"n\"><v>");
    push_g16(out, x);
    out.extend_from_slice(b"</v></c>");
}

/// Encode data rows 2..=n+1 of a sheet, calling `flush` with chunks of XML.
fn write_rows<F: FnMut(&[u8]) -> std::io::Result<()>>(cols: &[Column], mut flush: F) -> std::io::Result<()> {
    let n = cols.first().map_or(0, |c| c.len());
    let names: Vec<Vec<u8>> = (0..cols.len()).map(col_name).collect();
    let mut buf: Vec<u8> = Vec::with_capacity(1 << 20);
    let mut cref: Vec<u8> = Vec::with_capacity(16);
    for i in 0..n {
        let row = (i + 2).to_string();
        buf.extend_from_slice(b"<row r=\"");
        buf.extend_from_slice(row.as_bytes());
        buf.extend_from_slice(b"\">");
        for (j, col) in cols.iter().enumerate() {
            cref.clear();
            cref.extend_from_slice(&names[j]);
            cref.extend_from_slice(row.as_bytes());
            match col {
                Column::Int(v) => push_number(&mut buf, &cref, v[i] as f64),
                Column::Float(v) => {
                    let x = v[i];
                    if x.is_nan() {
                        push_inline_str(&mut buf, &cref, "");
                    } else if x.is_infinite() {
                        push_inline_str(&mut buf, &cref, if x > 0.0 { "inf" } else { "-inf" });
                    } else {
                        push_number(&mut buf, &cref, x);
                    }
                }
                Column::Bool(v) => {
                    buf.extend_from_slice(b"<c r=\"");
                    buf.extend_from_slice(&cref);
                    buf.extend_from_slice(if v[i] { b"\" t=\"b\"><v>1</v></c>" } else { b"\" t=\"b\"><v>0</v></c>" });
                }
                Column::Text(v) => push_inline_str(&mut buf, &cref, v[i].as_deref().unwrap_or("")),
            }
        }
        buf.extend_from_slice(b"</row>");
        if buf.len() > (1 << 20) - 4096 {
            flush(&buf)?;
            buf.clear();
        }
    }
    flush(&buf)
}

/// MS-DOS date and time of the current UTC time (zip entry timestamps).
fn dos_datetime() -> (u16, u16) {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0) as i64;
    let days = secs.div_euclid(86400);
    let rem = secs.rem_euclid(86400);
    // civil-from-days (Howard Hinnant)
    let z = days + 719468;
    let era = z.div_euclid(146097);
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    let date = (((y - 1980).clamp(0, 127) as u16) << 9) | ((m as u16) << 5) | d as u16;
    let time = (((rem / 3600) as u16) << 11) | ((((rem % 3600) / 60) as u16) << 5) | ((rem % 60) / 2) as u16;
    (time, date)
}

struct CentralEntry {
    name: String,
    crc: u32,
    csize: u64,
    usize: u64,
    offset: u64,
}

/// Write the workbook (a deflate-compressed zip archive of `parts`).
pub fn write_xlsx(path: &Path, parts: &[Part], sheets: &[Vec<Column>]) -> Result<(), String> {
    let io = |e: std::io::Error| format!("cannot write {}: {}", path.display(), e);
    let file = File::create(path).map_err(io)?;
    let mut w = BufWriter::with_capacity(1 << 22, file);
    let (time, date) = dos_datetime();
    let mut central: Vec<CentralEntry> = Vec::with_capacity(parts.len());
    for part in parts {
        let offset = w.stream_position().map_err(io)?;
        // local file header; crc and sizes are patched once the data is written
        let mut lh: Vec<u8> = Vec::with_capacity(30 + part.name.len());
        lh.extend_from_slice(&0x04034b50u32.to_le_bytes());
        lh.extend_from_slice(&20u16.to_le_bytes()); // version needed
        lh.extend_from_slice(&0u16.to_le_bytes()); // flags
        lh.extend_from_slice(&8u16.to_le_bytes()); // deflate
        lh.extend_from_slice(&time.to_le_bytes());
        lh.extend_from_slice(&date.to_le_bytes());
        lh.extend_from_slice(&[0u8; 12]); // crc, compressed size, uncompressed size
        lh.extend_from_slice(&(part.name.len() as u16).to_le_bytes());
        lh.extend_from_slice(&0u16.to_le_bytes()); // extra length
        lh.extend_from_slice(part.name.as_bytes());
        w.write_all(&lh).map_err(io)?;
        let data_start = w.stream_position().map_err(io)?;

        let mut crc = Crc::new();
        let mut usize_: u64 = 0;
        {
            let level = if part.sheet.is_some() { Compression::fast() } else { Compression::default() };
            let mut enc = DeflateEncoder::new(&mut w, level);
            let mut put = |bytes: &[u8]| -> std::io::Result<()> {
                crc.update(bytes);
                usize_ += bytes.len() as u64;
                enc.write_all(bytes)
            };
            put(&part.head).map_err(io)?;
            if let Some(s) = part.sheet {
                let cols = sheets.get(s).ok_or("sheet index out of range")?;
                write_rows(cols, &mut put).map_err(io)?;
            }
            put(&part.tail).map_err(io)?;
            enc.finish().map_err(io)?;
        }
        let end = w.stream_position().map_err(io)?;
        let csize = end - data_start;
        if csize > u32::MAX as u64 || usize_ > u32::MAX as u64 || offset > u32::MAX as u64 {
            return Err("workbook too large for the fast writer".into());
        }
        w.seek(SeekFrom::Start(offset + 14)).map_err(io)?;
        w.write_all(&crc.sum().to_le_bytes()).map_err(io)?;
        w.write_all(&(csize as u32).to_le_bytes()).map_err(io)?;
        w.write_all(&(usize_ as u32).to_le_bytes()).map_err(io)?;
        w.seek(SeekFrom::Start(end)).map_err(io)?;
        central.push(CentralEntry {
            name: part.name.clone(),
            crc: crc.sum(),
            csize,
            usize: usize_,
            offset,
        });
    }
    let cd_start = w.stream_position().map_err(io)?;
    for e in &central {
        let mut ch: Vec<u8> = Vec::with_capacity(46 + e.name.len());
        ch.extend_from_slice(&0x02014b50u32.to_le_bytes());
        ch.extend_from_slice(&20u16.to_le_bytes()); // version made by
        ch.extend_from_slice(&20u16.to_le_bytes()); // version needed
        ch.extend_from_slice(&0u16.to_le_bytes()); // flags
        ch.extend_from_slice(&8u16.to_le_bytes()); // deflate
        ch.extend_from_slice(&time.to_le_bytes());
        ch.extend_from_slice(&date.to_le_bytes());
        ch.extend_from_slice(&e.crc.to_le_bytes());
        ch.extend_from_slice(&(e.csize as u32).to_le_bytes());
        ch.extend_from_slice(&(e.usize as u32).to_le_bytes());
        ch.extend_from_slice(&(e.name.len() as u16).to_le_bytes());
        ch.extend_from_slice(&[0u8; 12]); // extra len, comment len, disk, internal attr, external attr (4)
        ch.extend_from_slice(&(e.offset as u32).to_le_bytes());
        ch.extend_from_slice(e.name.as_bytes());
        w.write_all(&ch).map_err(io)?;
    }
    let cd_end = w.stream_position().map_err(io)?;
    let mut eocd: Vec<u8> = Vec::with_capacity(22);
    eocd.extend_from_slice(&0x06054b50u32.to_le_bytes());
    eocd.extend_from_slice(&[0u8; 4]); // disk numbers
    eocd.extend_from_slice(&(central.len() as u16).to_le_bytes());
    eocd.extend_from_slice(&(central.len() as u16).to_le_bytes());
    eocd.extend_from_slice(&((cd_end - cd_start) as u32).to_le_bytes());
    eocd.extend_from_slice(&(cd_start as u32).to_le_bytes());
    eocd.extend_from_slice(&0u16.to_le_bytes()); // comment length
    w.write_all(&eocd).map_err(io)?;
    w.flush().map_err(io)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn g(x: f64) -> String {
        let mut v = Vec::new();
        push_g16(&mut v, x);
        String::from_utf8(v).unwrap()
    }

    #[test]
    fn percent_g16() {
        let cases: &[(f64, &str)] = &[
            (0.0, "0"),
            (-0.0, "-0"),
            (706.0, "706"),
            (1.0, "1"),
            (0.1, "0.1"),
            (2.4032636839886794, "2.403263683988679"),
            (0.054795863080029594, "0.05479586308002959"),
            (1234567890123456.5, "1234567890123456"),
            (1e16, "1e+16"),
            (10000000000000001.0, "1e+16"),
            (0.0001, "0.0001"),
            (0.00001, "1e-05"),
            (7.305745779742904e-181, "7.305745779742904e-181"),
            (-1.0, "-1"),
            (123456.0, "123456"),
            (1.5e300, "1.5e+300"),
        ];
        for &(x, want) in cases {
            assert_eq!(g(x), want, "{:e}", x);
        }
    }

    #[test]
    fn column_names() {
        assert_eq!(col_name(0), b"A");
        assert_eq!(col_name(25), b"Z");
        assert_eq!(col_name(26), b"AA");
        assert_eq!(col_name(27), b"AB");
        assert_eq!(col_name(701), b"ZZ");
        assert_eq!(col_name(702), b"AAA");
    }
}
