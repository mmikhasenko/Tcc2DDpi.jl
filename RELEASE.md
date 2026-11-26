# Release gates

1. Merge the resolution and X(3872) examples once they execute in their documented environments. The resolution demo plots threshold offsets in MeV and converts them with `e2m` for the mass-based phase-space API.
2. On the release commit, require passing Julia 1.12 package tests and a successful Pages build. Verify the generated index links to both new notebook exports and the README, and check the example outputs for execution errors.
3. Bump the package version and publish the release candidate from that verified commit. The mass-based API has package version `0.2.0`; the earlier `v0.1.1-rc1` target applied to the energy-API maintenance release.
4. Keep the mass-interface PR synchronized with the release baseline. Before marking it ready, migrate the example callers, test both legacy and mass-grid cache loading, and compare the nominal phase space and amplitude to the release baseline at the same physical masses. Record the units and tolerance used for the comparison.
