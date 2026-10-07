# Contributing to SiliconeLab

The project is starting with three real pilot workflows: silicone oil, PDMS
networks, and silicone coatings. See `workflows/pilots.json` for candidate cases
and `docs/pilot-acceptance.md` for validation expectations.

Use issues to report reproducible problems, discuss scientific requirements, or
propose interfaces. Include the relevant source version, configuration, command,
and a small input/output example when possible. Keep credentials and personal
cluster settings outside submitted examples.

For implementation changes, describe the scientific behavior and its validation
in a pull request. Preserve links to the source tools and make any changed
scientific assumptions explicit. The shared engine and per-system interfaces
will continue to be established through the pilots. The first imported component
is the generator; detailed configuration rules are in `docs/generator.md`.

Generator changes must pass the CMake/CTest checks in the README. Keep vendored
scientific files unchanged unless an intentional scientific change is documented.
When updating an upstream import, record its commit and imported-file checksums
in `upstream.json`, preserve its license, and rerun native parity checks.
