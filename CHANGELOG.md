# CHANGELOG

## 0.4.0 --- Epistatic Transaction (Jan 15, 2025)

- Changed where output is stored: everything is now stored in `--outdir`
- Paths contained in the settings file will now be loaded relative to the
  settings file
- Added epistatic interactions
- Added preliminary cached fitness evaluation for recombinants
- Huge performance improvements when collecting all mutations
- Tree will now merge internal nodes when possible
- Methods that access data have to be annotated with `#[require_deferred_drop]`
- Introduced generic encoding for symbolic values `virolution::encoding::Symbol`
- Use smallvec for storing mutations, this should reduce the number of unnecessary
  allocations especially when mutants go extinct right away
- Host trait: Support for custom host implementations
- Attribute API: Create custom attributes with lazy and eager evaluation

## 0.3.0 --- Mixed Performance (May 06, 2024)

- Supply diverse initial population
- Fixed, then introduced caching with multi host populations
- Fixed, removed unwanted existence of defective interfering particels
- In single threaded mode, the migration will now be computed more precisely,
  i.e. there should be no more rounding errors when during transfer and very
  small migration rates will cause sporadic migration depending on the size of
  their residue.

## 0.2.0 --- Treecherous Harmony (Sep 04, 2023)

- Change encoding order to ATCG
- Leaking the Historian
- Trees only report the `block_id`
- Include `block_id` in fasta output
- Truncate barcode file at the start of the simulation
- Change naming in configuration files for consistency
    - `SimulationSettings` -> `Settings`
    - `SimulationParameters` -> `Parameters`
    - `SimulationPlan` -> `Schedule`
- Field names have been renamed accordingly
- Moved all things related to `Settings` into config module
- Moved several structures into core module
- Move application control flow into runner module
- Added individual ancestry tracing (this is costly...)

## 0.1.0 --- Rapid Expansion (Jul 26, 2023)

- Start versioning virolution

