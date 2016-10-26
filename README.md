# gotoh
An implementation of the Gotoh alignment algorithm.

[![Build Status](https://travis-ci.org/isovic/gotoh.svg?branch=master)](https://travis-ci.org/isovic/gotoh)

** Project under construction! **  
Very unstable yet, use at your own risk.  
The only reason why the project is public is to mess around with TravisCI.  

The goal of this project is to create an alignment library which implements global, semiglobal and local alignment using Gotoh's approach.  
The library currently provides all these modes but in quadratic memory.  
Next step is to implement the Hirschberg algorithm to reduce memory consumption to linear.  
After that, SIMD support will be added.  
