# Install Linux Essentials
apt-get update
apt-get install build-essential ca-certificates -y
apt-get install vim git wget -y

# Install R
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
echo "deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/" > /etc/apt/sources.list.d/R.list
apt-get update && apt-get install r-base-dev -y

# Install Python
apt-get install python3 python3-pip python-is-python3 -y

# Install R/Python Linux dependencies
apt-get install libudunits2-dev libgdal-dev libcairo2-dev libmagick++-dev libfontconfig1-dev -y
