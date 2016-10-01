The Ruby script snapp_prep.rb prepares XML format input files for the software SNAPP (http://beast2.org/snapp/), given a phylip format SNP matrix, a table linking species IDs and specimen IDs, and a file specifying age constraints. Optionally, a starting tree can be provided (recommended) and the number of MCMC iterations used by SNAPP can be specified. Files prepared by this script use a particular combination of priors and operators optimized for SNAPP analyses that aim to estimate species divergence times (note that all population sizes are linked in these analyses). Detailed annotation explaining the choice of priors and operators will be included in the XML file (unless turned off with option '-n').

A manuscript on divergence-time estimation with SNAPP is in preparation.

To learn about the available options of this script, run with '-h' (or without any options):

`ruby snapp_prep.rb -h`

For questions and bug fixes, email to michaelmatschiner@mac.com
