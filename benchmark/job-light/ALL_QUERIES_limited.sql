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

IMPORT INTO cast_info DSV "benchmark/job-light/data/cast_info.csv" ROWS 10000000;
IMPORT INTO title DSV "benchmark/job-light/data/title.csv" ROWS 10000000; 
IMPORT INTO movie_companies DSV "benchmark/job-light/data/movie_companies_cleaned.csv" ROWS 10000000; 
IMPORT INTO movie_info_idx DSV "benchmark/job-light/data/movie_info_idx_cleaned.csv" ROWS 10000000;
IMPORT INTO movie_keyword DSV "benchmark/job-light/data/movie_keyword.csv" ROWS 10000000;
IMPORT INTO movie_info DSV "benchmark/job-light/data/movie_info.csv" ROWS 10000000; 

-- 1
SELECT COUNT(*) FROM movie_companies, title, movie_info_idx WHERE title.id=movie_companies.movie_id AND title.id=movie_info_idx.movie_id AND movie_info_idx.info_type_id=112 AND movie_companies.company_type_id=2;
-- 2
SELECT COUNT(*) FROM movie_companies, title, movie_info_idx WHERE title.id=movie_companies.movie_id AND title.id=movie_info_idx.movie_id AND movie_info_idx.info_type_id=113 AND movie_companies.company_type_id=2 AND title.production_year>2005 AND title.production_year<2010;
-- 3
SELECT COUNT(*) FROM movie_companies, title, movie_info_idx WHERE title.id=movie_companies.movie_id AND title.id=movie_info_idx.movie_id AND movie_info_idx.info_type_id=112 AND movie_companies.company_type_id=2 AND title.production_year>2010;
-- 4
SELECT COUNT(*) FROM movie_companies, title, movie_info_idx WHERE title.id=movie_companies.movie_id AND title.id=movie_info_idx.movie_id AND movie_info_idx.info_type_id=113 AND movie_companies.company_type_id=2 AND title.production_year>2000;
-- 5
SELECT COUNT(*) FROM movie_companies, title, movie_keyword WHERE title.id=movie_companies.movie_id AND title.id=movie_keyword.movie_id AND movie_keyword.keyword_id=117;
-- 6
SELECT COUNT(*) FROM title, movie_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.production_year>2005;
-- 7
SELECT COUNT(*) FROM title, movie_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.production_year>2010;
-- 8
SELECT COUNT(*) FROM title, movie_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.production_year>1990;
-- 9
SELECT COUNT(*) FROM title, movie_info_idx, movie_keyword WHERE title.id=movie_info_idx.movie_id AND title.id=movie_keyword.movie_id AND title.production_year>2005 AND movie_info_idx.info_type_id=101;
-- 10
SELECT COUNT(*) FROM title, movie_info_idx, movie_keyword WHERE title.id=movie_info_idx.movie_id AND title.id=movie_keyword.movie_id AND title.production_year>2010 AND movie_info_idx.info_type_id=101;
-- 11
SELECT COUNT(*) FROM title, movie_info_idx, movie_keyword WHERE title.id=movie_info_idx.movie_id AND title.id=movie_keyword.movie_id AND title.production_year>1990 AND movie_info_idx.info_type_id=101;
-- 12
SELECT COUNT(*) FROM title, movie_info, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.production_year>2005 AND movie_companies.company_type_id=2;
-- 13
SELECT COUNT(*) FROM title, movie_info, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.production_year>2010 AND movie_companies.company_type_id=2;
-- 14
SELECT COUNT(*) FROM title, movie_info, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.production_year>1990 AND movie_companies.company_type_id=2;
-- 15
SELECT COUNT(*) FROM movie_keyword, title, cast_info WHERE title.id=movie_keyword.movie_id AND title.id=cast_info.movie_id AND title.production_year>2010 AND movie_keyword.keyword_id=8200;
-- 16
SELECT COUNT(*) FROM movie_keyword, title, cast_info WHERE title.id=movie_keyword.movie_id AND title.id=cast_info.movie_id AND title.production_year>2014;
-- 17
SELECT COUNT(*) FROM movie_keyword, title, cast_info WHERE title.id=movie_keyword.movie_id AND title.id=cast_info.movie_id AND title.production_year>2014 AND movie_keyword.keyword_id=8200;
-- 18
SELECT COUNT(*) FROM movie_keyword, title, cast_info WHERE title.id=movie_keyword.movie_id AND title.id=cast_info.movie_id AND title.production_year>2000 AND movie_keyword.keyword_id=8200;
-- 19
SELECT COUNT(*) FROM movie_keyword, title, cast_info WHERE title.id=movie_keyword.movie_id AND title.id=cast_info.movie_id AND title.production_year>2000;
-- 20
SELECT COUNT(*) FROM cast_info, title WHERE title.id=cast_info.movie_id AND title.production_year>1980 AND title.production_year<1995;
-- 21
SELECT COUNT(*) FROM cast_info, title WHERE title.id=cast_info.movie_id AND title.production_year>1980 AND title.production_year<1984;
-- 22
SELECT COUNT(*) FROM cast_info, title WHERE title.id=cast_info.movie_id AND title.production_year>1980 AND title.production_year<2010;
-- 23
SELECT COUNT(*) FROM cast_info, title, movie_companies WHERE title.id=cast_info.movie_id AND title.id=movie_companies.movie_id AND cast_info.role_id=2;
-- 24
SELECT COUNT(*) FROM cast_info, title, movie_companies WHERE title.id=cast_info.movie_id AND title.id=movie_companies.movie_id AND cast_info.role_id=4;
-- 25
SELECT COUNT(*) FROM cast_info, title, movie_companies WHERE title.id=cast_info.movie_id AND title.id=movie_companies.movie_id AND cast_info.role_id=7;
-- 26
SELECT COUNT(*) FROM cast_info, title, movie_companies WHERE title.id=cast_info.movie_id AND title.id=movie_companies.movie_id AND title.production_year>2005 AND title.production_year<2015 AND cast_info.role_id=2;
-- 27
SELECT COUNT(*) FROM cast_info, title, movie_companies WHERE title.id=cast_info.movie_id AND title.id=movie_companies.movie_id AND title.production_year>2007 AND title.production_year<2010 AND cast_info.role_id=2;
-- 28
SELECT COUNT(*) FROM title, cast_info, movie_companies WHERE title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND title.production_year>2005 AND cast_info.role_id=1;
-- 29
SELECT COUNT(*) FROM title, cast_info, movie_companies WHERE title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND title.production_year>2010 AND cast_info.role_id=1;
-- 30
SELECT COUNT(*) FROM title, cast_info, movie_companies WHERE title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND title.production_year>1990;
-- 31
SELECT COUNT(*) FROM title, movie_keyword, movie_companies WHERE title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND movie_keyword.keyword_id=398 AND movie_companies.company_type_id=2 AND title.production_year>1950 AND title.production_year<2000;
-- 32
SELECT COUNT(*) FROM title, movie_keyword, movie_companies WHERE title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND movie_keyword.keyword_id=398 AND movie_companies.company_type_id=2;
-- 33
SELECT COUNT(*) FROM title, movie_keyword, movie_companies WHERE title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND title.production_year>1950;
-- 34
SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_info_idx.movie_id AND title.id=movie_companies.movie_id AND movie_info_idx.info_type_id=101 AND movie_info.info_type_id=3 AND title.production_year>2005 AND title.production_year<2008 AND movie_companies.company_type_id=2;
-- 35
SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_info_idx.movie_id AND title.id=movie_companies.movie_id AND movie_info_idx.info_type_id=113 AND movie_info.info_type_id=105;
-- 36
SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_info_idx.movie_id AND title.id=movie_companies.movie_id AND movie_info_idx.info_type_id=101 AND movie_info.info_type_id=3 AND title.production_year>2000 AND title.production_year<2010 AND movie_companies.company_type_id=2;
-- 37
SELECT COUNT(*) FROM title, movie_info, movie_companies, movie_info_idx WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.id=movie_info_idx.movie_id AND title.kind_id=1 AND movie_companies.company_type_id=2 AND movie_info_idx.info_type_id=101 AND movie_info.info_type_id=16;
-- 38
SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_info_idx.movie_id AND title.production_year>2010 AND title.kind_id=1 AND movie_info.info_type_id=8 AND movie_info_idx.info_type_id=101;
-- 39
SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_info_idx.movie_id AND title.kind_id=1 AND movie_info.info_type_id=8 AND movie_info_idx.info_type_id=101;
-- 40
SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_info_idx.movie_id AND title.production_year>2005 AND movie_info.info_type_id=8 AND movie_info_idx.info_type_id=101;
-- 41
SELECT COUNT(*) FROM title, movie_info, movie_companies, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND movie_info.info_type_id=16 AND title.production_year>2000;
-- 42
SELECT COUNT(*) FROM title, movie_info, movie_companies, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND movie_info.info_type_id=16 AND title.production_year>2005 AND title.production_year<2010;
-- 43
SELECT COUNT(*) FROM title, movie_info, movie_companies, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND movie_info.info_type_id=16 AND title.production_year>1990;
-- 44
SELECT COUNT(*) FROM cast_info, title, movie_keyword, movie_companies WHERE title.id=cast_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND movie_keyword.keyword_id=117;
-- 45
SELECT COUNT(*) FROM title, movie_info, movie_info_idx, cast_info WHERE title.id=movie_info.movie_id AND title.id=movie_info_idx.movie_id AND title.id=cast_info.movie_id AND movie_info.info_type_id=105 AND movie_info_idx.info_type_id=100;
-- 46
SELECT COUNT(*) FROM title, movie_info, movie_info_idx, cast_info WHERE title.id=movie_info.movie_id AND title.id=movie_info_idx.movie_id AND title.id=cast_info.movie_id AND movie_info.info_type_id=3 AND movie_info_idx.info_type_id=101 AND title.production_year>2008 AND title.production_year<2014;
-- 47
SELECT COUNT(*) FROM title, movie_info, movie_info_idx, cast_info WHERE title.id=movie_info.movie_id AND title.id=movie_info_idx.movie_id AND title.id=cast_info.movie_id AND movie_info.info_type_id=3 AND movie_info_idx.info_type_id=100;
-- 48
SELECT COUNT(*) FROM title, movie_info, movie_companies, cast_info WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND cast_info.role_id=2 AND movie_info.info_type_id=16 AND title.production_year>2005 AND title.production_year<2009;
-- 49
SELECT COUNT(*) FROM title, movie_info, movie_companies, cast_info WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND cast_info.role_id=2 AND movie_info.info_type_id=16;
-- 50
SELECT COUNT(*) FROM title, movie_info, movie_companies, cast_info WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND cast_info.role_id=2 AND movie_info.info_type_id=16 AND title.production_year>2000;
-- 51
SELECT COUNT(*) FROM title, cast_info, movie_keyword WHERE title.id=movie_keyword.movie_id AND title.id=cast_info.movie_id AND title.production_year>1950 AND title.kind_id=1;
-- 52
SELECT COUNT(*) FROM title, cast_info, movie_keyword WHERE title.id=movie_keyword.movie_id AND title.id=cast_info.movie_id AND title.production_year>2000 AND title.kind_id=1;
-- 53
SELECT COUNT(*) FROM title, movie_keyword, movie_companies, movie_info WHERE title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND title.id=movie_info.movie_id AND movie_keyword.keyword_id=398 AND movie_companies.company_type_id=2 AND title.production_year>1950 AND title.production_year<2000;
-- 54
SELECT COUNT(*) FROM title, movie_keyword, movie_companies, movie_info WHERE title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND title.id=movie_info.movie_id AND movie_keyword.keyword_id=398 AND movie_companies.company_type_id=2 AND title.production_year>2000 AND title.production_year<2010;
-- 55
SELECT COUNT(*) FROM title, movie_keyword, movie_companies, movie_info WHERE title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND title.id=movie_info.movie_id AND movie_keyword.keyword_id=398 AND movie_companies.company_type_id=2 AND title.production_year>1950 AND title.production_year<2010;
-- 56
SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_keyword, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_info_idx.movie_id AND title.id=movie_companies.movie_id AND title.production_year>2008 AND movie_info.info_type_id=8 AND movie_info_idx.info_type_id=101;
-- 57
SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_keyword, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_info_idx.movie_id AND title.id=movie_companies.movie_id AND title.production_year>2009 AND movie_info.info_type_id=8 AND movie_info_idx.info_type_id=101;
-- 58
SELECT COUNT(*) FROM title, movie_info, movie_companies, cast_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND title.id=movie_keyword.movie_id AND cast_info.role_id=2 AND movie_info.info_type_id=16 AND title.production_year>2010;
-- 59
SELECT COUNT(*) FROM title, movie_info, movie_companies, cast_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND title.id=movie_keyword.movie_id AND cast_info.role_id=2 AND movie_info.info_type_id=16 AND title.production_year>2010 AND movie_companies.company_id=22956;
-- 60
SELECT COUNT(*) FROM title, movie_info, movie_companies, cast_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND title.id=movie_keyword.movie_id AND cast_info.role_id=2 AND movie_info.info_type_id=16 AND title.production_year>2000;
-- 61
SELECT COUNT(*) FROM title, movie_info, movie_info_idx, cast_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_info_idx.movie_id AND title.id=cast_info.movie_id AND title.id=movie_keyword.movie_id AND movie_info.info_type_id=3 AND movie_info_idx.info_type_id=100;
-- 62
SELECT COUNT(*) FROM title, movie_info, movie_info_idx, cast_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_info_idx.movie_id AND title.id=cast_info.movie_id AND title.id=movie_keyword.movie_id AND movie_info.info_type_id=3 AND movie_info_idx.info_type_id=100 AND title.production_year>2010;
-- 63
SELECT COUNT(*) FROM title, cast_info, movie_keyword, movie_info_idx WHERE title.id=movie_keyword.movie_id AND title.id=cast_info.movie_id AND title.id=movie_info_idx.movie_id AND title.production_year>2000 AND title.kind_id=1 AND movie_info_idx.info_type_id=101;
-- 64
SELECT COUNT(*) FROM title, cast_info, movie_keyword, movie_info_idx WHERE title.id=movie_keyword.movie_id AND title.id=cast_info.movie_id AND title.id=movie_info_idx.movie_id AND title.production_year>2005 AND title.kind_id=1 AND movie_info_idx.info_type_id=101;
-- 65
SELECT COUNT(*) FROM title, movie_keyword, movie_companies, movie_info WHERE title.id=movie_keyword.movie_id AND title.id=movie_companies.movie_id AND title.id=movie_info.movie_id AND movie_keyword.keyword_id=398 AND movie_companies.company_type_id=2 AND title.production_year=1998;
-- 66
SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_keyword, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_info_idx.movie_id AND title.id=movie_companies.movie_id AND title.production_year>2000 AND movie_info.info_type_id=8 AND movie_info_idx.info_type_id=101;
-- 67
SELECT COUNT(*) FROM title, movie_info, movie_info_idx, movie_keyword, movie_companies WHERE title.id=movie_info.movie_id AND title.id=movie_keyword.movie_id AND title.id=movie_info_idx.movie_id AND title.id=movie_companies.movie_id AND title.production_year>2005 AND movie_info.info_type_id=8 AND movie_info_idx.info_type_id=101;
-- 68
SELECT COUNT(*) FROM title, movie_info, movie_companies, cast_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND title.id=movie_keyword.movie_id AND cast_info.role_id=2 AND movie_info.info_type_id=16 AND title.production_year>2000 AND title.production_year<2010 AND movie_keyword.keyword_id=7084;
-- 69
SELECT COUNT(*) FROM title, movie_info, movie_companies, cast_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_companies.movie_id AND title.id=cast_info.movie_id AND title.id=movie_keyword.movie_id AND cast_info.role_id=2 AND movie_info.info_type_id=16 AND title.production_year>2000 AND title.production_year<2005 AND movie_keyword.keyword_id=7084;
-- 70
SELECT COUNT(*) FROM title, movie_info, movie_info_idx, cast_info, movie_keyword WHERE title.id=movie_info.movie_id AND title.id=movie_info_idx.movie_id AND title.id=cast_info.movie_id AND title.id=movie_keyword.movie_id AND movie_info.info_type_id=3 AND movie_info_idx.info_type_id=100 AND title.production_year>2000;