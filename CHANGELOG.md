# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [2.0.0] - 2019-06-28

### Added

- New option to choose Mathematica FEM solver with `Method->"NDSolve"`.

### Changed

- Analysis time step is now adaptive (options `StartingStepSize` and `MaxStepSize`) instead of constant.
- Default FEM mesh topology is triangular.
- All FEM mesh topologies supported by `ElementMesh` are accepted.

### Removed

- Option `"NoTimeSteps"` is removed because time step length is now adaptive.

## [1.0.0] - 2019-03-07

### Changed

- Cosmetic code updates
  
## [0.3.1] - 2019-03-06

### Fixed

- Documentation tweaks to reduce its size
- Updated project layout to be more standard

## [0.3.0] - 2019-01-25

### Changed

- Expanded documentation

## [0.2.0] - 2018-12-05

### Added

- Very basic documentation
- Test support for multiple OS

## [0.1.0] - 2018-11-28

### Added

- Initial functionality

[Unreleased]: https://github.com/c3m-labs/HeatTrans/compare/v2.0.0...HEAD
[2.0.0]: https://github.com/c3m-labs/HeatTrans/compare/v1.0.0...v2.0.0
[1.0.0]: https://github.com/c3m-labs/HeatTrans/compare/v0.3.1...v1.0.0
[0.3.1]: https://github.com/c3m-labs/HeatTrans/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/c3m-labs/HeatTrans/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/c3m-labs/HeatTrans/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/c3m-labs/HeatTrans/releases/tag/v0.1.0
