//! Fast writer for `result_all.txt`, producing the same bytes as
//! `pandas.DataFrame.to_csv(path, sep="\t", index=False)` for the column types YACHT writes
//! (int64, float64, bool and text columns). The Python caller falls back to pandas for anything else.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// One output column.
pub enum Column<'a> {
    Int(&'a [i64]),
    Float(&'a [f64]),
    Bool(&'a [bool]),
    /// None = missing value (written as the empty string)
    Text(Vec<Option<String>>),
}

impl Column<'_> {
    pub fn len(&self) -> usize {
        match self {
            Column::Int(v) => v.len(),
            Column::Float(v) => v.len(),
            Column::Bool(v) => v.len(),
            Column::Text(v) => v.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Shortest round-trip decimal digits of a finite double: (negative, digits without leading/trailing zeros,
/// decimal exponent of the first digit). Uses Ryu, which (like Python's repr and numpy) breaks exact ties
/// between two shortest candidates towards the even digit.
fn shortest_digits(x: f64) -> (bool, Vec<u8>, i32) {
    let mut buf = ryu::Buffer::new();
    let s = buf.format_finite(x);
    let (neg, s) = match s.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, s),
    };
    let (mant, exp) = match s.split_once(['e', 'E']) {
        Some((m, e)) => (m, e.parse::<i32>().expect("bad exponent from ryu")),
        None => (s, 0),
    };
    let (int_part, frac_part) = mant.split_once('.').unwrap_or((mant, ""));
    let mut digits: Vec<u8> = int_part.bytes().chain(frac_part.bytes()).collect();
    // the decimal point sits after `point` digits of `digits`
    let mut point = int_part.len() as i32 + exp;
    let lead = digits.iter().take_while(|&&c| c == b'0').count();
    digits.drain(..lead);
    point -= lead as i32;
    while digits.last() == Some(&b'0') {
        digits.pop();
    }
    if digits.is_empty() {
        return (neg, vec![b'0'], 0);
    }
    (neg, digits, point - 1)
}

/// Append Python's `repr(x)` of a finite or infinite double (NaN must be handled by the caller).
/// Python (and numpy's `astype(str)`, which pandas uses) print the shortest digit string that round-trips,
/// in positional notation when the decimal exponent is in [-4, 16) and in scientific notation otherwise.
pub fn push_float_repr(out: &mut Vec<u8>, x: f64) {
    if x.is_infinite() {
        out.extend_from_slice(if x > 0.0 { b"inf" } else { b"-inf" });
        return;
    }
    let (neg, digits, exp) = shortest_digits(x);
    if neg {
        out.push(b'-');
    }
    let decpt = exp + 1; // position of the decimal point relative to the digit string
    if (-4..16).contains(&exp) {
        if decpt <= 0 {
            out.extend_from_slice(b"0.");
            out.extend(std::iter::repeat(b'0').take((-decpt) as usize));
            out.extend_from_slice(&digits);
        } else if decpt as usize >= digits.len() {
            out.extend_from_slice(&digits);
            out.extend(std::iter::repeat(b'0').take(decpt as usize - digits.len()));
            out.extend_from_slice(b".0");
        } else {
            out.extend_from_slice(&digits[..decpt as usize]);
            out.push(b'.');
            out.extend_from_slice(&digits[decpt as usize..]);
        }
    } else {
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
    }
}

/// Append a text field with csv.QUOTE_MINIMAL quoting for a tab-separated file (the caller guarantees that
/// the field contains no line breaks).
fn push_text(out: &mut Vec<u8>, s: &str) {
    if s.bytes().any(|c| c == b'\t' || c == b'"') {
        out.push(b'"');
        for c in s.bytes() {
            if c == b'"' {
                out.push(b'"');
            }
            out.push(c);
        }
        out.push(b'"');
    } else {
        out.extend_from_slice(s.as_bytes());
    }
}

/// True if `s` can be written by this writer (no line breaks, which the csv module treats specially).
pub fn text_supported(s: &str) -> bool {
    !s.bytes().any(|c| c == b'\n' || c == b'\r')
}

pub fn write_tsv(path: &Path, header: &[String], columns: &[Column], line_terminator: &str) -> Result<(), String> {
    if header.len() != columns.len() {
        return Err("header and columns differ in length".into());
    }
    let n = columns.first().map_or(0, |c| c.len());
    if columns.iter().any(|c| c.len() != n) {
        return Err("columns differ in length".into());
    }
    let file = File::create(path).map_err(|e| format!("cannot write {}: {}", path.display(), e))?;
    let mut w = BufWriter::with_capacity(1 << 20, file);
    let io = |e: std::io::Error| format!("cannot write {}: {}", path.display(), e);
    let mut line: Vec<u8> = Vec::with_capacity(1024);
    for (i, h) in header.iter().enumerate() {
        if i > 0 {
            line.push(b'\t');
        }
        push_text(&mut line, h);
    }
    line.extend_from_slice(line_terminator.as_bytes());
    w.write_all(&line).map_err(io)?;
    let mut ibuf = itoa_buffer();
    for row in 0..n {
        line.clear();
        for (i, col) in columns.iter().enumerate() {
            if i > 0 {
                line.push(b'\t');
            }
            match col {
                Column::Int(v) => line.extend_from_slice(format_int(&mut ibuf, v[row])),
                Column::Float(v) => {
                    let x = v[row];
                    if !x.is_nan() {
                        push_float_repr(&mut line, x);
                    }
                }
                Column::Bool(v) => line.extend_from_slice(if v[row] { b"True" } else { b"False" }),
                Column::Text(v) => {
                    if let Some(s) = &v[row] {
                        push_text(&mut line, s);
                    }
                }
            }
        }
        line.extend_from_slice(line_terminator.as_bytes());
        w.write_all(&line).map_err(io)?;
    }
    w.flush().map_err(io)?;
    Ok(())
}

fn itoa_buffer() -> [u8; 24] {
    [0u8; 24]
}

fn format_int(buf: &mut [u8; 24], v: i64) -> &[u8] {
    let mut i = buf.len();
    let neg = v < 0;
    let mut u = v.unsigned_abs();
    loop {
        i -= 1;
        buf[i] = b'0' + (u % 10) as u8;
        u /= 10;
        if u == 0 {
            break;
        }
    }
    if neg {
        i -= 1;
        buf[i] = b'-';
    }
    &buf[i..]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(x: f64) -> String {
        let mut v = Vec::new();
        push_float_repr(&mut v, x);
        String::from_utf8(v).unwrap()
    }

    #[test]
    fn python_repr_layout() {
        let cases: &[(f64, &str)] = &[
            (0.0, "0.0"),
            (-0.0, "-0.0"),
            (1.0, "1.0"),
            (706.0, "706.0"),
            (0.5, "0.5"),
            (0.0001, "0.0001"),
            (0.00001, "1e-05"),
            (1.5e-300, "1.5e-300"),
            (1e16, "1e+16"),
            (1e15, "1000000000000000.0"),
            (123456789012345678.0, "1.2345678901234568e+17"),
            (2.4032636839886794, "2.4032636839886794"),
            (7.305745779742904e-181, "7.305745779742904e-181"),
            (0.054795863080029594, "0.054795863080029594"),
            (-1.0, "-1.0"),
            (f64::INFINITY, "inf"),
            (f64::NEG_INFINITY, "-inf"),
            (5e-324, "5e-324"),
            (1.7976931348623157e308, "1.7976931348623157e+308"),
            (1743829569681555.2, "1743829569681555.2"),
            (133906845957463.12, "133906845957463.12"),
            (-34019466669101.062, "-34019466669101.062"),
            (0.1, "0.1"),
            (100.0, "100.0"),
            (1e22, "1e+22"),
            (1.5e-7, "1.5e-07"),
        ];
        for &(x, want) in cases {
            assert_eq!(r(x), want, "{:e}", x);
        }
    }

    #[test]
    fn ints() {
        let mut b = itoa_buffer();
        assert_eq!(format_int(&mut b, 0), b"0");
        assert_eq!(format_int(&mut b, -12), b"-12");
        assert_eq!(format_int(&mut b, i64::MIN), b"-9223372036854775808");
        assert_eq!(format_int(&mut b, i64::MAX), b"9223372036854775807");
    }
}
