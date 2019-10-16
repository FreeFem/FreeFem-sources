#!/bin/bash
 
# script to generate a xml report, make a visualization with the Test Results Analyser plugin 

grep 'global-test-result: XFAIL' examples/*/*trs  | sed -e "s/.trs::global-test-result/ /g ;  s/\//:/g ; s/ : /:/g " | cut -c10- >> result.txt
grep 'global-test-result: CPU' examples/*/*trs  | sed -e "s/.trs::global-test-result/ /g ;  s/\//:/g ; s/ : /:/g " | cut -c10- >> result.txt
grep 'global-test-result: SKIP' examples/*/*trs  | sed -e "s/.trs::global-test-result/ /g ;  s/\//:/g ; s/ : /:/g " | cut -c10- >> result.txt
grep 'global-test-result: PASS' examples/*/*trs  | sed -e "s/.trs::global-test-result/ /g ;  s/\//:/g ; s/ : /:/g " | cut -c10- >> result.txt

tests=0 failures=0 errors=0 skipped=0 xerrors=0 cpu=0
Gtests=0 Gfailures=0 Gerrors=0 Gskipped=0 Gxerrors=0 Gcpu=0
endclass=0
folder=toto
first=1

echo "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" >> report.xml
echo "<testsuites name=\"FreeFEM\" duration=\"0.000\">" >> report.xml
while IFS=: read classname name result 
do

if [ "$classname" != "$folder" ] && [ $first==1 ] ; then
if [ "$folder" != 'toto' ]; then echo "    </testsuite>" >> report.xml ; fi
sed -i -e "s/beginTest/    <testsuite name=\"$folder\" tests=\"$tests\" failures=\"$failures\" errors=\"$errors\" skipped=\"$skipped\" time=\"0.000\">/g" report.xml
echo "beginTest" >> report.xml
tests=0 failures=0 errors=0 skipped=0 xerrors=0
folder=$classname
first=1
else
first=0
fi

if [ "$result" == "PASS" ]; then
    tests=$((tests+1))
        Gtests=$((Gtests+1))
        echo "        <testcase name=\"$name\" classname=\"$classname\"/>" >> report.xml
elif [ "$result" == "SKIP" ]; then
     skipped=$((skipped+1))
         Gskipped=$((Gskipped+1))
         echo "        <testcase name=\"$name\" classname=\"$classname\">" >> report.xml
         echo "            <skipped> Skipped test, not library/plugin on the VM </skipped> " >> report.xml
         echo "        </testcase>" >> report.xml
 elif [ "$result" == "CPU" ]; then
     cpu=$((cpu+1))
         Gcpu=$((Gcpu+1))
         echo "        <testcase name=\"$name\" classname=\"$classname\">" >> report.xml
         echo "            <failure> Cputime limit exceeded </failure> " >> report.xml
         echo "        </testcase>" >> report.xml
elif [ "$result" == "XFAIL" ]; then
     errors=$((errors+1))
     Gerrors=$((Gerrors+1))
     echo "        <testcase name=\"$name\" classname=\"$classname\">" >> report.xml
     echo "                  <error> This is the end of edp script with the error:" >> report.xml
     cat -v examples/$classname/$name.log | tail -35 | sed -e "s/\&/\&#38;;/g" -e "s/</\&lt;/g" -e "s/>/\&gt;/g" -e "s/'/\&#39;/g" >> report.xml
     echo "                 </error> " >> report.xml
     echo "        </testcase>" >> report.xml
elif [ "$result" == "FAIL" ]; then
     failures=$((failures+1))
     Gfailures=$((Gfailures+1))
     echo "        <testcase name=\"$name\" classname=\"$classname\">" >> report.xml
     echo "                  <error> This is the end of edp script with the error:" >> report.xml
     cat -v examples/$classname/$name.log | tail -35 | sed -e "s/\&/\&#38;;/g" -e "s/</\&lt;/g" -e "s/>/\&gt;/g" -e "s/'/\&#39;/g" >> report.xml
     echo "                 </error> " >> report.xml
     echo "        </testcase>" >> report.xml
fi

done < result.txt 2>/dev/null

echo "    </testsuite>" >> report.xml
echo "</testsuites>" >> report.xml
sed -i -e "s/beginTest/    <testsuite name=\"$folder\" tests=\"$tests\" failures=\"$failures\" errors=\"$errors\" skipped=\"$skipped\" time=\"0.000\">/g" report.xml

sed -i -e "s/<testsuites name=\"FreeFEM\" duration=\"0.000\">/<testsuites name=\"FreeFEM\" tests=\"$Gtests\" failures=\"$Gfailures\" errors=\"$Gerrors\" skipped=\"$Gskipped\" duration=\"0.000\">/g" report.xml

mv report.xml etc/jenkins/resultForJenkins/report.xml
echo " ****************Make check results************ "
echo " PASS: $Gtests  SKIPPED: $Gskipped   XFAIL: $Gerrors  FAIL:$Gfailures  CPU: $Gcpu"

