use phf::phf_map;

pub static TRANSFERS: phf::Map<&'static str, &[[f64; 4]; 4]> = phf_map! {
    "migration_fwd" => &MIGRATION_FWD,
    "migration_rev" => &MIGRATION_REV,

};

const MIGRATION_FWD: [[f64; 4]; 4] = [
    [1., 0., 0., 0.],
    [0.2, 0.8, 0., 0.],
    [0., 0.2, 0.8, 0.],
    [0., 0., 0., 0.],
];
const MIGRATION_REV: [[f64; 4]; 4] = [
    [0.8, 0.2, 0., 0.],
    [0., 1., 0., 0.],
    [0., 0., 1., 0.],
    [0., 0., 0., 0.],
];
