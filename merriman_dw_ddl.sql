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
    bp int not null,
    A1 nvarchar(6),
    A2 nvarchar(6),
	primary key (variantID)
);

# dimParameter
# Records the parameter type (Odds ratio (OR), or beta-estimate (beta))
drop table if exists dimParameter;
create table dimParameter (
	parameterTypeID tinyint not null,
    parameterTypeDesc nvarchar(4),
    primary key (parameterTypeID)
);
insert into dimParameter
	select 1, 'beta'
    union
    select 2, 'OR';

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

# Experimental Conditions table
# Records the combination of case / control conditions
# which is unique to each experiment.
# NOTE:
#    currentExperimentId: in the context of UKBiobank, this is the unique combination of conditions
#    experimentID is simply a surrogate key that we are using for convenience in the data warehouse
drop table if exists dimExperiment;
create table dimExperiment (
    experimentID tinyint not null auto_increment,
    caseCondition varchar(20) not null,
    controlCondition varchar(20) not null,
    currentExperimentID int,
    datasetID int not null,
    primary key (experimentID),
    foreign key (datasetID) references dimSourceDataset(datasetID),
    constraint uniq_Experiment unique (datasetID, currentExperimentID)
);
insert into dimExperiment
    select (@Nexp := @Nexp + 1),
            myControls.Diuretic,
            myCases.All,
            (@CurExp := @CurExp + 1),
            1
    from (
            select 'Diuretic'
		    union select 'Not diuretic'
            union select 'Unknown diuretic'
	) myControls, (
            select 'All'
            union select 'Winnard'
            union select 'Hospital'
            union select 'Self-reported'
            union select 'Self + ULT'
            union select 'Male (All)'
    ) myCases,
    (select @Nexp := 0) nexp, (select @CurExp := 0) cexp

;

# Experiment-iteration table
# This table provides a mechanism to associate
# each cohort or controls with a specific experiment iteration
drop table if exists dimIteration;
create table dimIteration (
    experimentID tinyint not null,
    iterationID tinyint not null,
    primary key (experimentID, iterationID),
    foreign key (experimentID) references dimExperiment(experimentID)
);

# Person table
# Records information about the people / samples
# NOTE:
#    gender encoded as 1 / 2, 1= MALE, 2 = FEMALE (PLINK coding)
#    iid & fid are the cooresponding columns from PLINK
#    personID is a unique identifier (surrogate key) created within the database.
# This schema DOES NOT check the uniqueness of (iid, fid). This table is essentially a
# collection of unique FAM files, rather than unique people (i.e. the same person
# could be in more than one study, (more than 1 FAM file) and therefore, they will
# appear more than once in this table with DIFFERENT personIDs.
# The data import pipelines should be cognizant of this.
#drop table if exists dimPerson;
#create table dimPerson (
#    personID int not null auto_increment,
#    iid nvarchar(8),
#    fid varchar(8),
#    age int,
#    gender int,
#    ethnicity varchar(32),
#    primary key (personID)
#);
# Discussion with Murray - have remove the Age, gender, ethnicity columns from dimPerson
# This means that dimPerson is no longer necessary. Removing, and updating
# dimPeopleCase, dimPeopleControl as appropriate.


# GWAS Fact table
# NOTE:
#    majorAllele - indicates which allele is the major allele (saves a lot of headaches)
#        encode as {1, 2} corresponding to dimVariant{A1, A2}.
drop table if exists factGWAS;
create table factGWAS (
    iterationID tinyint not null,
    experimentID tinyint not null,
    variantID int not null,
    parameterType tinyint not null,
    parameterEstimate float,
    stdError float,
    pValue float,
    MAF_U float,
    MAF_A float,
    hwe_geno_U varchar(32),
    hwe_geno_A varchar(32),
    hwe_P_U float,
    hwe_P_A float,
    majorAllele tinyint not null,
    primary key (iterationID, experimentID, variantID),
    foreign key (experimentID, iterationID) references dimIteration(experimentID, iterationID),
    foreign key (variantID) references dimVariant(variantID),
    foreign key (parameterType) references dimParameter(parameterTypeID)
);

# Person - to - fact mapping tables
# PersonCase and PersonControl provide a mapping between the
# people and a particular iteration of an experiment.
drop table if exists dimPersonCase;
drop table if exists dimPersonControl;
# dimPersonCase maps between dimPerson and dimExperiment
create table dimPersonCase (
    experimentID tinyint not null,
    iid nvarchar(32) not null,
    fid nvarchar(32) not null,
    foreign key (experimentID) references dimExperiment(experimentID),
    primary key (experimentID),
    constraint uniq_case_groups unique (experimentID, iid, fid)
);
# dimPersonControl maps between dimPerson and dimIteration
create table dimPersonControl (
    iterationID tinyint not null,
    experimentID tinyint not null,
    iid nvarchar(32) not null,
    fid nvarchar(32) not null,
    foreign key (experimentID, iterationID) references dimIteration(experimentID, iterationID),
    primary key (experimentID, iterationID),
    constraint uniq_control_groups unique (experimentID, iterationID, iid, fid)
);


