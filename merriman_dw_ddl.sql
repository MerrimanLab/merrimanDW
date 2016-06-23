# DDL Script for merriman data warehouse
#
# Nick Burns
# April, 2016
drop database if exists merrimanDW;
create database merrimanDW;
use merrimanDW;

# Variants table
# Contains the Chromosome and Allele information for each SNP
# notes:  bp is the chromosomal position in bases
drop table if exists dimVariant;
create table dimVariant (
    variantID int auto_increment,
    chromosome tinyint not null,
    snp nvarchar(20) not null,
    GRCh37_bp int,
    GRCh38_bp int,
    A1 nvarchar(6),
    A2 nvarchar(6),
	primary key (variantID)
);

# dimSourceDataset
# Records meta data about the dataset being investigated
drop table if exists dimSourceDataset;
create table dimSourceDataset (
	datasetID int not null auto_increment,
    datasetName nvarchar(128) not null,
    platform nvarchar(128),
    genomeBuild nvarchar(8),
    primary key (datasetID)
);
insert into dimSourceDataset
select 1, 'UKBiobank', NULL, 'GRCh37';

# dimCondition
# CaseControl conditions e.g. gout, diabetic, obese... 
drop table if exists dimCondition;
create table dimCondition (
    conditionID int not null auto_increment,
    conditionDescription nvarchar(128) not null,
    primary key (conditionID)
);

# Experimental Conditions table
# Records meta-data specific to each GWAS experiment.
# NOTES:  
#    Decription: this is a TEXT field (i.e. free text field) which can be used 
#                to add additional information about the experiment
#    Covariates: also a free TEXT field, important to record the covariates used in the GWAS model.
#    ParameterType: 0 = 'odds ratio', 1 = 'regression (beta) estimate'
#    DataLoadDate: a DATETIME field, defaults to the DATETIME of insert.
drop table if exists dimExperiment;
create table dimExperiment (
    experimentID int not null auto_increment,
    trait nvarchar(32) not null,
    caseCondition int not null,
    controlCondition int not null,
    datasetID int not null,
    Description TEXT,
    Covariates TEXT not null,
    ParameterType int not null,
    DataLoadDate TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    primary key (experimentID),
    foreign key (datasetID) references dimSourceDataset(datasetID),
    foreign key (caseCondition) references dimCondition(conditionID),
    foreign key (controlCondition) references dimCondition(conditionID)
);


# Person table
# Records information about the people / samples
# This is essentially a .fam file (iid, fid, CaseControl)
# where CaseControl is: 0 = control, 1 = case
drop table if exists dimPerson;
create table dimPerson (
    personID int not null auto_increment,
    iid nvarchar(8),
    fid varchar(8),
    CaseControl tinyint not null,
    experimentID int not null,
    primary key (personID),
    foreign key (experimentID) references dimExperiment(experimentID)
);

# GWAS Fact table
# NOTE:
#    majorAllele - indicates which allele is the major allele (saves a lot of headaches)
#        encode as {1, 2} corresponding to dimVariant{A1, A2}.
drop table if exists factGWAS;
create table factGWAS (
    experimentID int not null,
    variantID int not null,
    parameterEstimate float,
    stdError float,
    pValue float,
    MAF_U float,
    MAF_A float,
    MAF_ALL float,
    hwe_geno_U varchar(32),
    hwe_geno_A varchar(32),
    hwe_geno_ALL varchar(32),
    hwe_P_U float,
    hwe_P_A float,
    hwe_P_ALL float,
    majorAllele tinyint not null,
    infoscore float,
    primary key (experimentID, variantID),
    foreign key (experimentID) references dimExperiment(experimentID),
    foreign key (variantID) references dimVariant(variantID)
);



