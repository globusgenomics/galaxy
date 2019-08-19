#!/bin/bash
version=1.4.3
mypath=`dirname $0`
svn update $mypath
rev=`svnversion $mypath`
if [ $rev == "exported" ]
	then
	echo "Error! The project does not seem to be under version control"
	exit 1
fi
rev=`echo $rev | sed 's/[A-Z]//g' | sed 's/[0-9]*://'`
echo "current revision is $rev"
rev=`expr $rev + 1`
echo "SVN status:"
svn status $mypath
read -p "continue with commit? (y/n) "
if [ "$REPLY" == "y" ]
	then
	echo "Setting version to $version and revision to $rev"
	sed -i -e "s/\\\$REVISION=[0-9]*/\$REVISION=$rev/" -e "s/\\\$VERSION=\"[0-9\.]*\"/\$VERSION=\"$version\"/" $mypath/Atlas-Indel2/Atlas-Indel2.rb $mypath/Atlas-SNP2/Atlas-SNP2.rb
	sed -i -e "s/\REVISION=\"[0-9]*\"/REVISION=\"$rev\"/" -e "s/VERSION=\"[0-9\.]*\"/VERSION=\"$version\"/" $mypath/vcfPrinter/settings.rb.default
	svn commit -m "$1" $mypath
fi
