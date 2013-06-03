# awk -f keys.awk examples++-tutorial/lestables|sort -u >ff-name
$2 > 256 && $2 < 1000 && NF == 3  && /^   / {  {print $1,"keywords"}}

/type :/ { if (nn !="") {print nn,"variable"}; nn= substr($1,1,length($1)-1); #print nn, "++";
 if(!nn) {nn= substr($2,1,length($2)-1);  } }
/operator/ && nn {op=substr($0,12,2); print nn,op ; nn ="" }  #   operator. : $
END {  if (nn !="") {print nn,"variable"} }
