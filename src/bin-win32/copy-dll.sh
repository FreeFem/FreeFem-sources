ldd *.exe |grep -v /WINDOWS/|grep -v /bin-win32/ |grep '=>' |awk '{print "cp",$3,"."}'|sort -u|grep -v '[?]'|sh
