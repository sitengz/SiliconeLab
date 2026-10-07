# Operating scope

Cluster access, remote directory writes, submission, monitoring, downloads,
recovery, and extra sampling are distinct capabilities. User authorization and
the cluster's access policy establish which are available for a particular run.

Record approved directories, resource limits, concurrency, and recovery rules in
the user profile and run specification. Once established, routine actions within
that scope should not trigger repeated confirmation.

Changing a force field, temperature, composition, reaction target, analysis
interval, or simulation duration changes the scientific specification. Such
changes require either an explicit pre-established rule or the researcher's
decision, and must be recorded.

Authentication requirements are unresolved until the user's connection procedure
is inspected. Never assume SSH keys, VPN access, MFA behavior, or unattended-login
permission from the presence of a cluster hostname.
