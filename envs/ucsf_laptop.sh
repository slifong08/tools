# sarahfong
 
# updated
# 20190109 vgi01 alias
# 20190403 scpdat function added
 
 
# description
# bash profile for sarah's desktop

alias python='python3'
 
# ACCRE
alias accre='ssh login.accre.vanderbilt.edu'
alias capra1='ssh capra1.accre.vanderbilt.edu'
 
function nb1 { ssh -L 8888:capra1:"$1" fongsl@capra1.accre.vanderbilt.edu; }
 
export -f nb1
 
# WYNTON
 
function wnb1 { ssh fongsl@dev2.wynton.ucsf.edu -J fongsl@log1.wynton.ucsf.edu -L "$2":localhost:"$1"; }

export -f wnb1

function wnbr { ssh fongsl@devr8.wynton.ucsf.edu -J fongsl@log1.wynton.ucsf.edu -L "$2":localhost:"$1"; }

export -f wnbr

function gnb1 { ssh fongsl@gpudev1.wynton.ucsf.edu -J fongsl@log1.wynton.ucsf.edu -L "$2":localhost:"$1"; }

export -f gnb1

function gnb8 { ssh fongsl@gpudevr8.wynton.ucsf.edu -J fongsl@log1.wynton.ucsf.edu -L "$2":localhost:"$1"; }

export -f gnb1
 
 

# virtual environments

# python 3.6, ete trees, clustalo
#alias phylo='conda activate /Users/sarahfong/miniconda3/envs/phylo'

# basic pandas, scipy, sklearn, etc. 
alias venv='micromamba activate /Users/sarahfong/micromamba/envs/mambaenv'

#export PATH

#PYTOOLS = "$HOME/tools/py_"

export PYTOOLS

USR = "/usr/lib"
export USR


eval "$(/opt/homebrew/bin/brew shellenv)"
alias wynton="ssh dev2.wynton.ucsf.edu"
alias gpu="ssh gpudev1.wynton.ucsf.edu"

alias log1="ssh log1.wynton.ucsf.edu"
alias log2="ssh log2.wynton.ucsf.edu"

alias mamba="micromamba"
