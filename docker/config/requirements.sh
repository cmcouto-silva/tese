# Plink
wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20201019.zip
unzip plink_linux_x86_64_20201019.zip -d ./plink
mv plink/plink /usr/bin/
rm -rf plink*

# Plink2
wget http://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20210123.zip
unzip plink2_linux_x86_64_20210123.zip
mv plink2 /usr/bin/
rm -rf plink2*

# Admixture
wget http://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz
tar xvzf admixture_linux-1.3.0.tar.gz
cp dist/admixture_linux-1.3.0/admixture /usr/bin/
rm admixture_linux-1.3.0.tar.gz
