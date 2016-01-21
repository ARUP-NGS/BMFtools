set -e
cd dmp && python hashdmp_test.py && cd ..
cd marksplit && python marksplit_test.py && cd ..
