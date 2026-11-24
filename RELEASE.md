# Release gates

1. Merge the resolution and X(3872) examples once they execute in their documented environments. The resolution demo uses threshold offsets in MeV, as does the current interpolated-channel API.
2. On the release commit, require passing Julia 1.12 package tests and a successful Pages build. Verify the generated index links to both new notebook exports and the README, and check the example outputs for execution errors.
3. Bump the package version and publish the release candidate from that verified commit. Use `v0.1.1-rc1` for the existing API maintenance release.
4. Rebase the mass-interface PR onto the release baseline. Before marking it ready, migrate the example callers, test both legacy and mass-grid cache loading, and compare the nominal phase space and amplitude to the release baseline at the same physical masses. Record the units and tolerance used for the comparison.
