use std::io::{Error, Write};

use serde::Serialize;

#[derive(Serialize)]
pub struct BarcodeEntry<'a> {
    pub barcode: &'a String,
    pub experiment: &'a String,
    pub time: usize,
    pub replicate: usize,
    pub compartment: usize,
}

impl BarcodeEntry<'_> {
    // pub fn new(
    //     barcode: &String,
    //     experiment: &String,
    //     time: usize,
    //     replicate: usize,
    //     compartment: usize,
    // ) -> Self {
    //     Self {
    //         barcode,
    //         experiment,
    //         time,
    //         replicate,
    //         compartment,
    //     }
    // }

    pub fn write(&self, writer: &mut impl Write) -> Result<(), Error> {
        let mut csv_writer = csv::Writer::from_writer(writer);
        csv_writer.serialize(self)?;
        Ok(())
    }
}
