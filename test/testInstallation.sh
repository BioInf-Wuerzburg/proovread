#!/bin/bash

testInstallation() {
    ../bin/proovread --pre test --overwrite --long ../sample/F.antasticus_long_error.fq --short ../sample/F.antasticus_short.fq --sr-qv-offset 33 &>/dev/null
    returnCode=$?
    assertEquals "Testrun exited with exit code of 0" 0 $returnCode
    assertTrue "test.trimmed.fq exists" "[ -e test/test.trimmed.fq ]"
    assertTrue "test.trimmed.fa exists" "[ -e test/test.trimmed.fa ]"
    assertTrue "test.untrimmed.fq exists" "[ -e test/test.untrimmed.fq ]"
}

. /tmp/shunit2/shunit2
