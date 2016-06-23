# DDL script for selection database
#
# Murray Cadzow
# June 2016
drop database if exists selectionDW;
create database selectionDW;
use selectionDW;

# position table
# Contains the Chromosome and Allele information for each SNP
# notes:  bp is the chromosomal position in bases
drop table if exists dimPos;
create table dimPos (
    posID int not null auto_increment,
    chrom int,
    chrom_start int,
    chrom_end int,
    slide int,
	primary key (posID)
);

# dimPopData
# contains the information about the populations
drop table if exists dimPopData;
create table dimPopData (
	popID int not null auto_increment,
    pop nvarchar(10) not null,
    description nvarchar(140),
    superPop nvarchar(3),
    genotypePlatform varchar(32),
	primary key (popID)
);

# dimStat
# selection statistic names 
drop table if exists dimStat;
create table dimStat (
    statID int not null auto_increment,
    statName nvarchar(20) not null,
    statDescription nvarchar(128),
    primary key (statID)
);

# dimGene
# positions of genes by ensembl ID or gene name 
drop table if exists dimGene;
create table dimGene (
    posID int not null,
    GeneName nvarchar(32) not null,
    EnsGeneID nvarchar(32),
    primary key (EnsGeneID, GeneName),
    foreign key (posID) references dimPos(posID)
);

# dimExperiment
drop table if exists dimExperiment;
create table dimExperiment (
	experimentID int not null auto_increment,
    description varchar(200),
    primary key (experimentID)
);



# intraSel
# store the inter population stats
drop table if exists intraSel;
create table intraSel (
    posID int not null,
    popID int not null,
    statValue float,
    statID int not null,
    experimentID int not null,
    primary key (popID, posID, statID),
    foreign key (posID) references dimPos(posID),
    foreign key (popID) references dimPopData(popID),
    foreign key (statID) references dimStat(statID),
    foreign key (experimentID) references dimExperiment(experimentID)
);

# intraSel
# store the inter population stats
drop table if exists interSel;
create table interSel (
    posID int not null,
    popID1 int not null,
    popID2 int not null,
    statValue float,
    statID int not null,
    experimentID int not null,
    foreign key (posID) references dimPos(posID),
    foreign key (popID1) references dimPopData(popID),
    foreign key (popID2) references dimPopData(popID),
    foreign key (statID) references dimStat(statID),
    foreign key (experimentID) references dimExperiment(experimentID)
);
