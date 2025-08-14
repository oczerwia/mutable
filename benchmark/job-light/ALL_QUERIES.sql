CREATE DATABASE job_light;
USE job_light;

CREATE TABLE cast_info (
    id INT(4) NOT NULL PRIMARY KEY,
    person_id INT(4) NOT NULL,
    movie_id INT(4) NOT NULL,
    person_role_id INT(4),
    note CHAR(100),
    nr_order INT(4),
    role_id INT(4) NOT NULL
);

CREATE TABLE title (
    id INT(4) NOT NULL PRIMARY KEY,
    title CHAR(100) NOT NULL,
    imdb_index CHAR(5),
    kind_id INT(4) NOT NULL,
    production_year INT(4),
    imdb_id INT(4),
    phonetic_code CHAR(5),
    episode_of_id INT(4),
    season_nr INT(4),
    episode_nr INT(4),
    series_years CHAR(49),
    md5sum CHAR(32)
);



CREATE TABLE movie_companies (
    id INT(4) NOT NULL PRIMARY KEY,
    movie_id INT(4) NOT NULL,
    company_id INT(4) NOT NULL,
    company_type_id INT(4) NOT NULL,
    note CHAR(128)
);


CREATE TABLE movie_info_idx (
    id INT(4) NOT NULL PRIMARY KEY,
    movie_id INT(4) NOT NULL,
    info_type_id INT(4) NOT NULL,
    info CHAR(100) NOT NULL,
    note CHAR(128)
);

CREATE TABLE movie_keyword (
    id INT(4) NOT NULL PRIMARY KEY,
    movie_id INT(4) NOT NULL,
    keyword_id INT(4) NOT NULL
);

CREATE TABLE movie_info (
    id INT(4) NOT NULL PRIMARY KEY,
    movie_id INT(4) NOT NULL,
    info_type_id INT(4) NOT NULL,
    info CHAR(128) NOT NULL,
    note CHAR(128)
);



-- Import data from CSV files
IMPORT INTO cast_info DSV "benchmark/job-light/data/cast_info.csv" ROWS 10000;
IMPORT INTO title DSV "benchmark/job-light/data/title.csv" ROWS 10000;
IMPORT INTO movie_companies DSV "benchmark/job-light/data/movie_companies_cleaned.csv" ROWS 10000;
IMPORT INTO movie_info_idx DSV "benchmark/job-light/data/movie_info_idx_cleaned.csv" ROWS 10000;
IMPORT INTO movie_keyword DSV "benchmark/job-light/data/movie_keyword.csv" ROWS 10000;
IMPORT INTO movie_info DSV "benchmark/job-light/data/movie_info.csv" ROWS 10000;




SELECT COUNT(*) FROM movie_companies, title, movie_info_idx WHERE title.id=movie_companies.movie_id AND title.id=movie_info_idx.movie_id AND movie_info_idx.info_type_id=112 AND movie_companies.company_type_id=2;

SELECT COUNT(*) FROM movie_companies, title, movie_info_idx WHERE title.id=movie_companies.movie_id AND title.id=movie_info_idx.movie_id AND movie_info_idx.info_type_id=113 AND movie_companies.company_type_id=2 AND title.production_year>2005 AND title.production_year<2010;

SELECT COUNT(*) FROM movie_companies, title, movie_info_idx WHERE title.id=movie_companies.movie_id AND title.id=movie_info_idx.movie_id AND movie_info_idx.info_type_id=112 AND movie_companies.company_type_id=2 AND title.production_year>2010;

SELECT COUNT(*) FROM movie_companies, title, movie_info_idx WHERE title.id=movie_companies.movie_id AND title.id=movie_info_idx.movie_id AND movie_info_idx.info_type_id=113 AND movie_companies.company_type_id=2 AND title.production_year>2000;

SELECT COUNT(*) FROM movie_companies, title, movie_keyword WHERE title.id=movie_companies.movie_id AND title.id=movie_keyword.movie_id AND movie_keyword.keyword_id=117;

SELECT COUNT(*) FROM title, movie_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.production_year>2005;

SELECT COUNT(*) FROM title, movie_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.production_year>2010;

SELECT COUNT(*) FROM title, movie_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.production_year>1990;

SELECT COUNT(*) FROM title, movie_info_idx, movie_keyword WHERE title.id=movie_info_idx.movie_id AND title.id=movie_keyword.movie_id AND title.production_year>2005 AND movie_info_idx.info_type_id=101;

SELECT COUNT(*) FROM title, movie_info_idx, movie_keyword WHERE title.id=movie_info_idx.movie_id AND title.id=movie_keyword.movie_id AND title.production_year>2010 AND movie_info_idx.info_type_id=101;

SELECT COUNT(*) FROM title, movie_info_idx, movie_keyword WHERE title.id=movie_info_idx.movie_id AND title.id=movie_keyword.movie_id AND title.production_year>1990 AND movie_info_idx.info_type_id=101;

SELECT COUNT(*) FROM title, movie_info, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.production_year>2005 AND movie_companies.company_type_id=2;

SELECT COUNT(*) FROM title, movie_info, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.production_year>2010 AND movie_companies.company_type_id=2;

SELECT COUNT(*) FROM title, movie_info, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.production_year>1990 AND movie_companies.company_type_id=2;

SELECT COUNT(*) FROM movie_keyword, title, cast_info WHERE title.id=movie_keyword.movie_id AND title.id=cast_info.movie_id AND title.production_year>2010 AND movie_keyword.keyword_id=8200;

SELECT COUNT(*) FROM movie_keyword, title, cast_info WHERE title.id=movie_keyword.movie_id AND title.id=cast_info.movie_id AND title.production_year>2014;

SELECT COUNT(*) FROM movie_keyword, title, cast_info WHERE title.id=movie_keyword.movie_id AND title.id=cast_info.movie_id AND title.production_year>2014 AND movie_keyword.keyword_id=8200;

SELECT COUNT(*) FROM movie_keyword, title, cast_info WHERE title.id=movie_keyword.movie_id AND title.id=cast_info.movie_id AND title.production_year>2000 AND movie_keyword.keyword_id=8200;

SELECT COUNT(*) FROM movie_keyword, title, cast_info WHERE title.id=movie_keyword.movie_id AND title.id=cast_info.movie_id AND title.production_year>2000;

SELECT COUNT(*) FROM cast_info, title WHERE title.id=cast_info.movie_id AND title.production_year>1980 AND title.production_year<1995;

SELECT COUNT(*) FROM cast_info, title WHERE title.id=cast_info.movie_id AND title.production_year>1980 AND title.production_year<1984;

SELECT COUNT(*) FROM cast_info, title WHERE title.id=cast_info.movie_id AND title.production_year>1980 AND title.production_year<2010;

SELECT COUNT(*) FROM cast_info, title, movie_companies WHERE title.id=cast_info.movie_id AND title.id=movie_companies.movie_id AND cast_info.role_id=2;

SELECT COUNT(*) FROM cast_info, title, movie_companies WHERE title.id=cast_info.movie_id AND title.id=movie_companies.movie_id AND cast_info.role_id=4;

SELECT COUNT(*) FROM cast_info, title, movie_companies WHERE title.id=cast_info.movie_id AND title.id=movie_companies.movie_id AND cast_info.role_id=7;

SELECT COUNT(*) FROM cast_info, title, movie_companies WHERE title.id=cast_info.movie_id AND title.id=movie_companies.movie_id AND title.production_year>2005 AND title.production_year<2015 AND cast_info.role_id=2;

SELECT COUNT(*) FROM cast_info, title, movie_companies WHERE title.id=cast_info.movie_id AND title.id=movie_companies.movie_id AND title.production_year>2007 AND title.production_year<2010 AND cast_info.role_id=2;

SELECT COUNT(*) FROM title, cast_info, movie_companies WHERE title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND title.production_year>2005 AND cast_info.role_id=1;

SELECT COUNT(*) FROM title, cast_info, movie_companies WHERE title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND title.production_year>2010 AND cast_info.role_id=1;

SELECT COUNT(*) FROM title, cast_info, movie_companies WHERE title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND title.production_year>1990;

SELECT COUNT(*) FROM title, movie_keyword, movie_companies WHERE title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND movie_keyword.keyword_id=398 AND movie_companies.company_type_id=2 AND title.production_year>1950 AND title.production_year<2000;

SELECT COUNT(*) FROM title, movie_keyword, movie_companies WHERE title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND movie_keyword.keyword_id=398 AND movie_companies.company_type_id=2;

SELECT COUNT(*) FROM title, movie_keyword, movie_companies WHERE title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND title.production_year>1950;

SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_info_idx.movie_id AND title.id=movie_companies.movie_id AND movie_info_idx.info_type_id=101 AND movie_info.info_type_id=3 AND title.production_year>2005 AND title.production_year<2008 AND movie_companies.company_type_id=2;

SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_info_idx.movie_id AND title.id=movie_companies.movie_id AND movie_info_idx.info_type_id=113 AND movie_info.info_type_id=105;

SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_info_idx.movie_id AND title.id=movie_companies.movie_id AND movie_info_idx.info_type_id=101 AND movie_info.info_type_id=3 AND title.production_year>2000 AND title.production_year<2010 AND movie_companies.company_type_id=2;

SELECT COUNT(*) FROM title, movie_info, movie_companies, movie_info_idx WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.id=movie_info_idx.movie_id AND title.kind_id=1 AND movie_companies.company_type_id=2 AND movie_info_idx.info_type_id=101 AND movie_info.info_type_id=16;

SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_info_idx.movie_id AND title.production_year>2010 AND title.kind_id=1 AND movie_info.info_type_id=8 AND movie_info_idx.info_type_id=101;

SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_info_idx.movie_id AND title.kind_id=1 AND movie_info.info_type_id=8 AND movie_info_idx.info_type_id=101;

SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_info_idx.movie_id AND title.production_year>2005 AND movie_info.info_type_id=8 AND movie_info_idx.info_type_id=101;

SELECT COUNT(*) FROM title, movie_info, movie_companies, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND movie_info.info_type_id=16 AND title.production_year>2000;

SELECT COUNT(*) FROM title, movie_info, movie_companies, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND movie_info.info_type_id=16 AND title.production_year>2005 AND title.production_year<2010;

SELECT COUNT(*) FROM title, movie_info, movie_companies, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND movie_info.info_type_id=16 AND title.production_year>1990;

SELECT COUNT(*) FROM cast_info, title, movie_keyword, movie_companies WHERE title.id=cast_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND movie_keyword.keyword_id=117;

SELECT COUNT(*) FROM title, movie_info, movie_info_idx, cast_info WHERE title.id=movie_info.movie_id AND title.id=movie_info_idx.movie_id AND title.id=cast_info.movie_id AND movie_info.info_type_id=105 AND movie_info_idx.info_type_id=100;

SELECT COUNT(*) FROM title, movie_info, movie_info_idx, cast_info WHERE title.id=movie_info.movie_id AND title.id=movie_info_idx.movie_id AND title.id=cast_info.movie_id AND movie_info.info_type_id=3 AND movie_info_idx.info_type_id=101 AND title.production_year>2008 AND title.production_year<2014;

SELECT COUNT(*) FROM title, movie_info, movie_info_idx, cast_info WHERE title.id=movie_info.movie_id AND title.id=movie_info_idx.movie_id AND title.id=cast_info.movie_id AND movie_info.info_type_id=3 AND movie_info_idx.info_type_id=100;

SELECT COUNT(*) FROM title, movie_info, movie_companies, cast_info WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND cast_info.role_id=2 AND movie_info.info_type_id=16 AND title.production_year>2005 AND title.production_year<2009;

SELECT COUNT(*) FROM title, movie_info, movie_companies, cast_info WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND cast_info.role_id=2 AND movie_info.info_type_id=16;

SELECT COUNT(*) FROM title, movie_info, movie_companies, cast_info WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND cast_info.role_id=2 AND movie_info.info_type_id=16 AND title.production_year>2000;

SELECT COUNT(*) FROM title, cast_info, movie_keyword WHERE title.id=movie_keyword.movie_id AND title.id=cast_info.movie_id AND title.production_year>1950 AND title.kind_id=1;

SELECT COUNT(*) FROM title, cast_info, movie_keyword WHERE title.id=movie_keyword.movie_id AND title.id=cast_info.movie_id AND title.production_year>2000 AND title.kind_id=1;

SELECT COUNT(*) FROM title, movie_keyword, movie_companies, movie_info WHERE title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND title.id=movie_info.movie_id AND movie_keyword.keyword_id=398 AND movie_companies.company_type_id=2 AND title.production_year>1950 AND title.production_year<2000;

SELECT COUNT(*) FROM title, movie_keyword, movie_companies, movie_info WHERE title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND title.id=movie_info.movie_id AND movie_keyword.keyword_id=398 AND movie_companies.company_type_id=2 AND title.production_year>2000 AND title.production_year<2010;

SELECT COUNT(*) FROM title, movie_keyword, movie_companies, movie_info WHERE title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND title.id=movie_info.movie_id AND movie_keyword.keyword_id=398 AND movie_companies.company_type_id=2 AND title.production_year>1950 AND title.production_year<2010;

SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_keyword, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_info_idx.movie_id AND title.id=movie_companies.movie_id AND title.production_year>2008 AND movie_info.info_type_id=8 AND movie_info_idx.info_type_id=101;

SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_keyword, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_info_idx.movie_id AND title.id=movie_companies.movie_id AND title.production_year>2009 AND movie_info.info_type_id=8 AND movie_info_idx.info_type_id=101;

SELECT COUNT(*) FROM title, movie_info, movie_companies, cast_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND title.id=movie_keyword.movie_id AND cast_info.role_id=2 AND movie_info.info_type_id=16 AND title.production_year>2010;

SELECT COUNT(*) FROM title, movie_info, movie_companies, cast_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND title.id=movie_keyword.movie_id AND cast_info.role_id=2 AND movie_info.info_type_id=16 AND title.production_year>2010 AND movie_companies.company_id=22956;

SELECT COUNT(*) FROM title, movie_info, movie_companies, cast_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND title.id=movie_keyword.movie_id AND cast_info.role_id=2 AND movie_info.info_type_id=16 AND title.production_year>2000;

SELECT COUNT(*) FROM title, movie_info, movie_info_idx, cast_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_info_idx.movie_id AND title.id=cast_info.movie_id AND title.id=movie_keyword.movie_id AND movie_info.info_type_id=3 AND movie_info_idx.info_type_id=100;

SELECT COUNT(*) FROM title, movie_info, movie_info_idx, cast_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_info_idx.movie_id AND title.id=cast_info.movie_id AND title.id=movie_keyword.movie_id AND movie_info.info_type_id=3 AND movie_info_idx.info_type_id=100 AND title.production_year>2010;

SELECT COUNT(*) FROM title, cast_info, movie_keyword, movie_info_idx WHERE title.id=movie_keyword.movie_id AND title.id=cast_info.movie_id AND title.id=movie_info_idx.movie_id AND title.production_year>2000 AND title.kind_id=1 AND movie_info_idx.info_type_id=101;

SELECT COUNT(*) FROM title, cast_info, movie_keyword, movie_info_idx WHERE title.id=movie_keyword.movie_id AND title.id=cast_info.movie_id AND title.id=movie_info_idx.movie_id AND title.production_year>2005 AND title.kind_id=1 AND movie_info_idx.info_type_id=101;

SELECT COUNT(*) FROM title, movie_keyword, movie_companies, movie_info WHERE title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND title.id=movie_info.movie_id AND movie_keyword.keyword_id=398 AND movie_companies.company_type_id=2 AND title.production_year=1998;

SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_keyword, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_info_idx.movie_id AND title.id=movie_companies.movie_id AND title.production_year>2000 AND movie_info.info_type_id=8 AND movie_info_idx.info_type_id=101;

SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_keyword, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_info_idx.movie_id AND title.id=movie_companies.movie_id AND title.production_year>2005 AND movie_info.info_type_id=8 AND movie_info_idx.info_type_id=101;

SELECT COUNT(*) FROM title, movie_info, movie_companies, cast_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND title.id=movie_keyword.movie_id AND cast_info.role_id=2 AND movie_info.info_type_id=16 AND title.production_year>2000 AND title.production_year<2010 AND movie_keyword.keyword_id=7084;

SELECT COUNT(*) FROM title, movie_info, movie_companies, cast_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND title.id=movie_keyword.movie_id AND cast_info.role_id=2 AND movie_info.info_type_id=16 AND title.production_year>2000 AND title.production_year<2005 AND movie_keyword.keyword_id=7084;

SELECT COUNT(*) FROM title, movie_info, movie_info_idx, cast_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_info_idx.movie_id AND title.id=cast_info.movie_id AND title.id=movie_keyword.movie_id AND movie_info.info_type_id=3 AND movie_info_idx.info_type_id=100 AND title.production_year>2000;