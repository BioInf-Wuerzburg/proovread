#!/bin/bash

testMultipleSeqs() {
    ../bin/proovread --pre test_multiple_seqs --overwrite --long ../sample/F.antasticus_long_error.fq --short ../sample/F.antasticus_short.fq --sr-qv-offset 33 2>log_multiple_seqs
    returnCode=$?
    assertEquals "Testrun exited with exit code of 0" 0 $returnCode
    numWarnings=$(grep -c "Use of uninitialized value \$bpN in division" log_multiple_seqs)
    assertEquals "No warnings messages are given indicating a problem with SeqFilter" 0 $numWarnings
}

testSingleSeq() {
    ../bin/proovread --pre test_single_seq --overwrite --long ../sample/F.antasticus_long_error.singleseq.fq --short ../sample/F.antasticus_short.fq --sr-qv-offset 33 2>log_single_seq
    returnCode=$?
    assertEquals "Testrun exited with exit code of 0" 0 $returnCode
    numWarnings=$(grep -c "Use of uninitialized value \$bpN in division" log_single_seq)
    assertEquals "No warnings messages are given indicating a problem with SeqFilter" 0 $numWarnings
}

. /tmp/shunit2/shunit2
