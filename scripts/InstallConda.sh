
# Links ####
# installing Conda: https://engineeringfordatascience.com/posts/install_miniconda_from_the_command_line/
# installing feems: https://github.com/NovembreLab/feems


# Installing Miniconda ####

mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh


# Installing feems using the Bioconda recipe ####

conda install -c bioconda feems -c conda-forge


# Installing pandas_plink

conda install -c conda-forge pandas-plink

# in the requirements of feems: pandas-plink==2.0.4
# But, I did not manage to install it with the following command because of pacakge dependencies: conda install -c conda-forge pandas-plink==2.0.4 

# Requirement networkx==2.4.0
conda install -c conda-forge networkx==2.4.0


# Running the python script with the feems analyses:
python scripts/feeems.py # being in the H2Pinpin directory

