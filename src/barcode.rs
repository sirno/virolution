use std::io::{Error, Write};

use serde::Serialize;

// TODO: Move to readwrite module and add documentation
#[derive(Serialize)]
pub struct BarcodeEntry<'a> {
    pub barcode: &'a String,
    pub experiment: &'a String,
    pub time: usize,
    pub replicate: usize,
    pub compartment: usize,
}

impl BarcodeEntry<'_> {
    pub fn write(&self, writer: &mut impl Write) -> Result<(), Error> {
        let mut csv_writer = csv::WriterBuilder::new()
            .has_headers(false)
            .from_writer(writer);
        csv_writer.serialize(self)?;
        csv_writer.flush()?;
        Ok(())
    }

    pub fn write_header(writer: &mut impl Write) -> Result<(), Error> {
        let mut csv_writer = csv::WriterBuilder::new()
            .has_headers(false)
            .from_writer(writer);
        csv_writer.write_record(["barcode", "experiment", "time", "replicate", "compartment"])?;
        csv_writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write() {
        let mut buffer = Vec::new();
        let entry = BarcodeEntry {
            barcode: &String::from("barcode"),
            experiment: &String::from("experiment"),
            time: 0,
            replicate: 0,
            compartment: 0,
        };
        BarcodeEntry::write_header(&mut buffer).unwrap();
        entry.write(&mut buffer).unwrap();
        entry.write(&mut buffer).unwrap();
        assert_eq!(
            buffer,
            b"barcode,experiment,time,replicate,compartment\nbarcode,experiment,0,0,0\nbarcode,experiment,0,0,0\n"
        )
    }
}
