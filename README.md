The code for the first problem is in P1/src/main.rs; the other files in P1 are just there to ensure the project builds correctly on machines other than my own. (It should build with any recent version of Rust. After installing Rust, execute "cargo run" inside the P1 directory.)

P2.m requires the MATLAB Parallel Computing Toolbox, but you can just delete the lines "parpool;" and "delete(gcp);" and change the "parfor" keyword to "for" to make it run on vanilla MATLAB. The simulation is pretty slow and parallelization is straightforward though.
