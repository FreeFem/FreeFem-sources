# Installation of specific packages - Not distributed
# ---------------------------------------------------

# This file is not distributed because it contains sensitive
# information about ftp servers.

# $Id$

# WWW access
# ----------

www=rascasse.inria.fr:www_gamma/Gamma/cdrom/ftp/freefem/
ftp=rascasse.inria.fr:/ftp_gamma/freefem/
ftpk=baobab.ann.jussieu.fr:public_html/ftp/freefem/

www:  FreeFem++v$(VERSION)_Win.zip  FreeFem++v$(VERSION)_MacOS.sit FreeFem++v$(VERSION)_MacOsX.tgz
	scp  HISTORY DOC/manual.pdf DOC/manual.ps.gz freefem++-$(VERSION).tar.gz FreeFem++v$(VERSION)_MacOsX.tgz FreeFem++v$(VERSION)_Win.zip  FreeFem++v$(VERSION)_MacOS.sit  $(ftpk)
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf FreeFem++v$(VERSION)_Win.zip  freefem++.zip"
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf FreeFem++v$(VERSION)_MacOS.sit freefem++.sit"
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf freefem++-$(VERSION).tar.gz freefem++.tar.gz"
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf freefem++-$(VERSION).tar.gz freefem++.tgz"
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf FreeFem++v$(VERSION)_MacOsX.tgz freefem++_MacOsX.tgz"
	scp  HISTORY  baobab.ann.jussieu.fr:www/.
	ssh  baobab.ann.jussieu.fr "cd www/.; ./.update++ $(VERSION) <ff++.htmx >freefem++.htm"
	scp  HISTORY DOC/manual.pdf DOC/manual.ps.gz freefem++-$(VERSION).tar.gz FreeFem++v$(VERSION)_MacOsX.tgz FreeFem++v$(VERSION)_Win.zip  FreeFem++v$(VERSION)_MacOS.sit  $(ftp)
	scp HISTORY DOC/manual.pdf DOC/manual.ps.gz  $(www)
	scp FreeFem++v$(VERSION)_Win.zip $(www)/freefem++.zip
	scp freefem++-$(VERSION).tar.gz $(www)/freefem++.tar.gz
	scp FreeFem++v$(VERSION)_MacOS.sit $(www)/freefem++.sit
	ssh rascasse.inria.fr "cd /ftp_gamma/freefem;ln -sf FreeFem++v$(VERSION)_Win.zip  freefem++.zip"
	ssh rascasse.inria.fr "cd /ftp_gamma/freefem;ln -sf FreeFem++v$(VERSION)_MacOS.sit freefem++.sit"
	ssh rascasse.inria.fr "cd /ftp_gamma/freefem;ln -sf freefem++-$(VERSION).tar.gz freefem++.tar.gz"
	ssh rascasse.inria.fr "cd /ftp_gamma/freefem;ln -sf freefem++-$(VERSION).tar.gz freefem++.tgz"
	ssh rascasse.inria.fr "cd /ftp_gamma/freefem;ln -sf FreeFem++v$(VERSION)_MacOsX.tgz freefem++_MacOsX.tgz"
	ssh rascasse.inria.fr "cd ~/public_html/; ./.update++ $(VERSION) <ff++.htmx >freefem++.htm"
