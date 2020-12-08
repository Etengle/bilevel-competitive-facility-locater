#!/bin/bash
./BCFL -i --o=o.txt | tee i.log
./BCFL --i=o.txt | tee o.log
diff i.log o.log | tee cmp.log
diff cmp.log io_cmp.golden
if [ "$?" == 0 ]; then
    echo "PASS!!!"
else
    echo "FAILED!!!"
fi
