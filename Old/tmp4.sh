cat 04.launch_whites_topmed.sh | sed "s/transformation/$1/" | sed "s/protein/$2/" | sed "s/visit/$3/" > ../results/${1}${2}_Visit${3}_W/launch_whites_topmed.sh
bash ../results/${1}${2}_Visit${3}_W/launch_whites_topmed.sh
#sleep 2h

