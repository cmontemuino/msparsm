[![Stories in Ready](https://badge.waffle.io/cmontemuino/mspar.png?label=in%20progress&title=In%20Progress)](https://waffle.io/cmontemuino/mspar)
mspar [![Build Status](https://travis-ci.org/cmontemuino/mspar.svg?branch=master)](https://travis-ci.org/cmontemuino/mspar)
=====

Parallel version of "ms" coalescent simulator using a master worker approach and a MPI implementation with on-demand scheduling.

# Test
In the **tests/cases** folder there is a set of test cases that can be used for performance testing.

Running the test case defined by *params.case.01*: `sh testcase.sh 01`

Output is going to be put in a new folder named as case**xx**, where the *xx* corresponds to the case number. Following the previous example,
the new folder is going to be named as **case01**.

The *testcase.sh* script is going to run 7 times. If you need less, then just updated the script. It is adviced to take a look over the script in
order to know its restrictions.
