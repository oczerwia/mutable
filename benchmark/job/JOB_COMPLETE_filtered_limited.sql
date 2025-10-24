CREATE DATABASE job;
USE job;

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


----- NEW


CREATE TABLE aka_name (
    id INT(4) NOT NULL PRIMARY KEY,
    person_id INT(4) NOT NULL,
    name CHAR(100),
    imdb_index CHAR(3),
    name_pcode_cf CHAR(11),
    name_pcode_nf CHAR(11),
    surname_pcode CHAR(11),
    md5sum CHAR(65)
);

CREATE TABLE aka_title (
    id INT(4) NOT NULL PRIMARY KEY,
    movie_id INT(4) NOT NULL,
    title CHAR(100),
    imdb_index CHAR(4),
    kind_id INT(4) NOT NULL,
    production_year INT(4),
    phonetic_code CHAR(5),
    episode_of_id INT(4),
    season_nr INT(4),
    episode_nr INT(4),
    note CHAR(72),
    md5sum CHAR(32)
);

CREATE TABLE char_name (
    id INT(4) NOT NULL PRIMARY KEY,
    name CHAR(100) NOT NULL,
    imdb_index CHAR(2),
    imdb_id INT(4),
    name_pcode_nf CHAR(5),
    surname_pcode CHAR(5),
    md5sum CHAR(32)
);

CREATE TABLE comp_cast_type (
    id INT(4) NOT NULL PRIMARY KEY,
    kind CHAR(32) NOT NULL
);

CREATE TABLE company_name (
    id INT(4) NOT NULL PRIMARY KEY,
    name CHAR(100) NOT NULL,
    country_code CHAR(6),
    imdb_id INT(4),
    name_pcode_nf CHAR(5),
    name_pcode_sf CHAR(5),
    md5sum CHAR(32)
);

CREATE TABLE company_type (
    id INT(4) NOT NULL PRIMARY KEY,
    kind CHAR(32)
);

CREATE TABLE complete_cast (
    id INT(4) NOT NULL PRIMARY KEY,
    movie_id INT(4),
    subject_id INT(4) NOT NULL,
    status_id INT(4) NOT NULL
);

CREATE TABLE info_type (
    id INT(4) NOT NULL PRIMARY KEY,
    info CHAR(32) NOT NULL
);

CREATE TABLE keyword (
    id INT(4) NOT NULL PRIMARY KEY,
    keyword CHAR(100) NOT NULL,
    phonetic_code CHAR(5)
);

CREATE TABLE kind_type (
    id INT(4) NOT NULL PRIMARY KEY,
    kind CHAR(15)
);

CREATE TABLE link_type (
    id INT(4) NOT NULL PRIMARY KEY,
    link CHAR(32) NOT NULL
);

CREATE TABLE movie_link (
    id INT(4) NOT NULL PRIMARY KEY,
    movie_id INT(4) NOT NULL,
    linked_movie_id INT(4) NOT NULL,
    link_type_id INT(4) NOT NULL
);

CREATE TABLE name (
    id INT(4) NOT NULL PRIMARY KEY,
    name CHAR(100) NOT NULL,
    imdb_index CHAR(9),
    imdb_id INT(4),
    gender CHAR(1),
    name_pcode_cf CHAR(5),
    name_pcode_nf CHAR(5),
    surname_pcode CHAR(5),
    md5sum CHAR(32)
);

CREATE TABLE role_type (
    id INT(4) NOT NULL PRIMARY KEY,
    role CHAR(32) NOT NULL
);

CREATE TABLE person_info (
    id INT(4) NOT NULL PRIMARY KEY,
    person_id INT(4) NOT NULL,
    info_type_id INT(4) NOT NULL,
    info CHAR(128) NOT NULL,
    note CHAR(128)
);

IMPORT INTO cast_info DSV "benchmark/job-data/cast_info.csv" ROWS 1000;
IMPORT INTO title DSV "benchmark/job-data/title.csv" ROWS 1000; 
IMPORT INTO movie_companies DSV "benchmark/job-data/movie_companies_cleaned.csv" ROWS 1000; 
IMPORT INTO movie_info_idx DSV "benchmark/job-data/movie_info_idx_cleaned.csv" ROWS 1000;
IMPORT INTO movie_keyword DSV "benchmark/job-data/movie_keyword.csv" ROWS 1000;
IMPORT INTO movie_info DSV "benchmark/job-data/movie_info.csv" ROWS 1000; 
IMPORT INTO aka_name DSV "benchmark/job-data/aka_name.csv" ROWS 1000;
IMPORT INTO aka_title DSV "benchmark/job-data/aka_title.csv" ROWS 1000;
IMPORT INTO char_name DSV "benchmark/job-data/char_name.csv" ROWS 1000;
IMPORT INTO comp_cast_type DSV "benchmark/job-data/comp_cast_type.csv" ROWS 1000;
IMPORT INTO company_name DSV "benchmark/job-data/company_name.csv" ROWS 1000;
IMPORT INTO company_type DSV "benchmark/job-data/company_type.csv" ROWS 1000;
IMPORT INTO complete_cast DSV "benchmark/job-data/complete_cast.csv" ROWS 1000;
IMPORT INTO info_type DSV "benchmark/job-data/info_type.csv" ROWS 1000;
IMPORT INTO keyword DSV "benchmark/job-data/keyword.csv" ROWS 1000;
IMPORT INTO kind_type DSV "benchmark/job-data/kind_type.csv" ROWS 1000;
IMPORT INTO link_type DSV "benchmark/job-data/link_type.csv" ROWS 1000;
IMPORT INTO movie_link DSV "benchmark/job-data/movie_link.csv" ROWS 1000;
IMPORT INTO name DSV "benchmark/job-data/name.csv" ROWS 1000;
IMPORT INTO role_type DSV "benchmark/job-data/role_type.csv" ROWS 1000;
IMPORT INTO person_info DSV "benchmark/job-data/person_info.csv" ROWS 1000;

-- Query 1
-- 10a
SELECT char_name.name, title.title FROM char_name, cast_info, company_name, company_type, movie_companies, role_type, title WHERE cast_info.note LIKE "%(voice)%" AND cast_info.note LIKE "%(uncredited)%" AND company_name.country_code = "[ru]" AND role_type.role = "actor" AND title.production_year > 2005 AND title.id = movie_companies.movie_id AND title.id = cast_info.movie_id AND cast_info.movie_id = movie_companies.movie_id AND char_name.id = cast_info.person_role_id AND role_type.id = cast_info.role_id AND company_name.id = movie_companies.company_id AND company_type.id = movie_companies.company_type_id;

-- Query 2
-- 10b
SELECT char_name.name, title.title FROM char_name, cast_info, company_name, company_type, movie_companies, role_type, title WHERE cast_info.note LIKE "%(producer)%" AND company_name.country_code = "[ru]" AND role_type.role = "actor" AND title.production_year > 2010 AND title.id = movie_companies.movie_id AND title.id = cast_info.movie_id AND cast_info.movie_id = movie_companies.movie_id AND char_name.id = cast_info.person_role_id AND role_type.id = cast_info.role_id AND company_name.id = movie_companies.company_id AND company_type.id = movie_companies.company_type_id;

-- Query 3
-- 10c
SELECT char_name.name, title.title FROM char_name, cast_info, company_name, company_type, movie_companies, role_type, title WHERE cast_info.note LIKE "%(producer)%" AND company_name.country_code = "[us]" AND title.production_year > 1990 AND title.id = movie_companies.movie_id AND title.id = cast_info.movie_id AND cast_info.movie_id = movie_companies.movie_id AND char_name.id = cast_info.person_role_id AND role_type.id = cast_info.role_id AND company_name.id = movie_companies.company_id AND company_type.id = movie_companies.company_type_id;

-- Query 4
-- 11a
SELECT company_name.name, link_type.link, title.title FROM company_name, company_type, keyword, link_type, movie_companies, movie_keyword, movie_link, title WHERE (company_name.name LIKE "%Film%" OR company_name.name LIKE "%Warner%") AND company_type.kind ="production companies" AND keyword.keyword ="sequel" AND link_type.link LIKE "%follow%" AND title.production_year >= 1950 AND title.production_year <= 2000 AND link_type.id = movie_link.link_type_id AND movie_link.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_type_id = company_type.id AND movie_companies.company_id = company_name.id AND movie_link.movie_id = movie_keyword.movie_id AND movie_link.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = movie_companies.movie_id;

-- Query 5
-- 11b
SELECT company_name.name, link_type.link, title.title FROM company_name, company_type, keyword, link_type, movie_companies, movie_keyword, movie_link, title WHERE (company_name.name LIKE "%Film%" OR company_name.name LIKE "%Warner%") AND company_type.kind ="production companies" AND keyword.keyword ="sequel" AND link_type.link LIKE "%follows%" AND title.production_year = 1998 AND title.title LIKE "%Money%" AND link_type.id = movie_link.link_type_id AND movie_link.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_type_id = company_type.id AND movie_companies.company_id = company_name.id AND movie_link.movie_id = movie_keyword.movie_id AND movie_link.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = movie_companies.movie_id;

-- Query 6
-- 11c
SELECT company_name.name, movie_companies.note, title.title FROM company_name, company_type, keyword, link_type, movie_companies, movie_keyword, movie_link, title WHERE (company_name.name LIKE "20th Century Fox%" OR company_name.name LIKE "Twentieth Century Fox%") AND (keyword.keyword = "sequel" OR keyword.keyword = "revenge" OR keyword.keyword = "based-on-novel") AND title.production_year > 1950 AND link_type.id = movie_link.link_type_id AND movie_link.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_type_id = company_type.id AND movie_companies.company_id = company_name.id AND movie_link.movie_id = movie_keyword.movie_id AND movie_link.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = movie_companies.movie_id;

-- Query 7
-- 11d
SELECT company_name.name, movie_companies.note, title.title FROM company_name, company_type, keyword, link_type, movie_companies, movie_keyword, movie_link, title WHERE (keyword.keyword = "sequel" OR keyword.keyword = "revenge" OR keyword.keyword = "based-on-novel") AND title.production_year > 1950 AND link_type.id = movie_link.link_type_id AND movie_link.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_type_id = company_type.id AND movie_companies.company_id = company_name.id AND movie_link.movie_id = movie_keyword.movie_id AND movie_link.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = movie_companies.movie_id;

-- Query 8
-- 12a
SELECT company_name.name, movie_info_idx.info, title.title FROM company_name, company_type, info_type, info_type, movie_companies, movie_info, movie_info_idx, title WHERE company_name.country_code = "[us]" AND company_type.kind = "production companies" AND info_type.info = "genres" AND info_type.info = "rating" AND (movie_info.info = "Drama" OR movie_info.info = "Horror") AND movie_info_idx.info > "8.0" AND title.production_year >= 2005 AND title.production_year <= 2008 AND title.id = movie_info.movie_id AND title.id = movie_info_idx.movie_id AND movie_info.info_type_id = info_type.id AND movie_info_idx.info_type_id = info_type.id AND title.id = movie_companies.movie_id AND company_type.id = movie_companies.company_type_id AND company_name.id = movie_companies.company_id AND movie_companies.movie_id = movie_info.movie_id AND movie_companies.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_info_idx.movie_id;

-- Query 9
-- 12b
SELECT movie_info.info, title.title FROM company_name, company_type, info_type, info_type, movie_companies, movie_info, movie_info_idx, title WHERE company_name.country_code ="[us]" AND (company_type.kind ="production companies" OR company_type.kind = "distributors") AND info_type.info ="budget" AND info_type.info ="bottom 10 rank" AND title.production_year >2000 AND (title.title LIKE "Birdemic%" OR title.title LIKE "%Movie%") AND title.id = movie_info.movie_id AND title.id = movie_info_idx.movie_id AND movie_info.info_type_id = info_type.id AND movie_info_idx.info_type_id = info_type.id AND title.id = movie_companies.movie_id AND company_type.id = movie_companies.company_type_id AND company_name.id = movie_companies.company_id AND movie_companies.movie_id = movie_info.movie_id AND movie_companies.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_info_idx.movie_id;

-- Query 10
-- 12c
SELECT company_name.name, movie_info_idx.info, title.title FROM company_name, company_type, info_type, info_type, movie_companies, movie_info, movie_info_idx, title WHERE company_name.country_code = "[us]" AND company_type.kind = "production companies" AND info_type.info = "genres" AND info_type.info = "rating" AND (movie_info.info = "Drama" OR movie_info.info = "Horror" OR movie_info.info = "Western" OR movie_info.info = "Family") AND movie_info_idx.info > "7.0" AND title.production_year >= 2000 AND title.production_year <= 2010 AND title.id = movie_info.movie_id AND title.id = movie_info_idx.movie_id AND movie_info.info_type_id = info_type.id AND movie_info_idx.info_type_id = info_type.id AND title.id = movie_companies.movie_id AND company_type.id = movie_companies.company_type_id AND company_name.id = movie_companies.company_id AND movie_companies.movie_id = movie_info.movie_id AND movie_companies.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_info_idx.movie_id;

-- Query 11
-- 13a
SELECT movie_info.info, movie_info_idx.info, title.title FROM company_name, company_type, info_type, info_type, kind_type, movie_companies, movie_info, movie_info_idx, title WHERE company_name.country_code ="[de]" AND company_type.kind ="production companies" AND info_type.info ="rating" AND info_type.info ="release dates" AND kind_type.kind ="movie" AND movie_info.movie_id = title.id AND info_type.id = movie_info.info_type_id AND kind_type.id = title.kind_id AND movie_companies.movie_id = title.id AND company_name.id = movie_companies.company_id AND company_type.id = movie_companies.company_type_id AND movie_info_idx.movie_id = title.id AND info_type.id = movie_info_idx.info_type_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_info_idx.movie_id = movie_companies.movie_id;

-- Query 12
-- 13b
SELECT company_name.name, movie_info_idx.info, title.title FROM company_name, company_type, info_type, info_type, kind_type, movie_companies, movie_info, movie_info_idx, title WHERE company_name.country_code ="[us]" AND company_type.kind ="production companies" AND info_type.info ="rating" AND info_type.info ="release dates" AND kind_type.kind ="movie" AND (title.title LIKE "%Champion%" OR title.title LIKE "%Loser%") AND movie_info.movie_id = title.id AND info_type.id = movie_info.info_type_id AND kind_type.id = title.kind_id AND movie_companies.movie_id = title.id AND company_name.id = movie_companies.company_id AND company_type.id = movie_companies.company_type_id AND movie_info_idx.movie_id = title.id AND info_type.id = movie_info_idx.info_type_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_info_idx.movie_id = movie_companies.movie_id;

-- Query 13
-- 13c
SELECT company_name.name, movie_info_idx.info, title.title FROM company_name, company_type, info_type, info_type, kind_type, movie_companies, movie_info, movie_info_idx, title WHERE company_name.country_code ="[us]" AND company_type.kind ="production companies" AND info_type.info ="rating" AND info_type.info ="release dates" AND kind_type.kind ="movie" AND (title.title LIKE "Champion%" OR title.title LIKE "Loser%") AND movie_info.movie_id = title.id AND info_type.id = movie_info.info_type_id AND kind_type.id = title.kind_id AND movie_companies.movie_id = title.id AND company_name.id = movie_companies.company_id AND company_type.id = movie_companies.company_type_id AND movie_info_idx.movie_id = title.id AND info_type.id = movie_info_idx.info_type_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_info_idx.movie_id = movie_companies.movie_id;

-- Query 14
-- 13d
SELECT company_name.name, movie_info_idx.info, title.title FROM company_name, company_type, info_type, info_type, kind_type, movie_companies, movie_info, movie_info_idx, title WHERE company_name.country_code ="[us]" AND company_type.kind ="production companies" AND info_type.info ="rating" AND info_type.info ="release dates" AND kind_type.kind ="movie" AND movie_info.movie_id = title.id AND info_type.id = movie_info.info_type_id AND kind_type.id = title.kind_id AND movie_companies.movie_id = title.id AND company_name.id = movie_companies.company_id AND company_type.id = movie_companies.company_type_id AND movie_info_idx.movie_id = title.id AND info_type.id = movie_info_idx.info_type_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_info_idx.movie_id = movie_companies.movie_id;

-- Query 15
-- 14a
SELECT movie_info_idx.info, title.title FROM info_type, info_type, keyword, kind_type, movie_info, movie_info_idx, movie_keyword, title WHERE info_type.info = "countries" AND info_type.info = "rating" AND (keyword.keyword = "murder" OR keyword.keyword = "murder-in-title" OR keyword.keyword = "blood" OR keyword.keyword = "violence") AND kind_type.kind = "movie" AND (movie_info.info = "Sweden" OR movie_info.info = "Norway" OR movie_info.info = "Germany" OR movie_info.info = "Denmark" OR movie_info.info = "Swedish" OR movie_info.info = "Denish" OR movie_info.info = "Norwegian" OR movie_info.info = "German" OR movie_info.info = "USA" OR movie_info.info = "American") AND movie_info_idx.info < "8.5" AND title.production_year > 2010 AND kind_type.id = title.kind_id AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_info_idx.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id;

-- Query 16
-- 14b
SELECT movie_info_idx.info, title.title FROM info_type, info_type, keyword, kind_type, movie_info, movie_info_idx, movie_keyword, title WHERE info_type.info = "countries" AND info_type.info = "rating" AND (keyword.keyword = "murder" OR keyword.keyword = "murder-in-title") AND kind_type.kind = "movie" AND (movie_info.info = "Sweden" OR movie_info.info = "Norway" OR movie_info.info = "Germany" OR movie_info.info = "Denmark" OR movie_info.info = "Swedish" OR movie_info.info = "Denish" OR movie_info.info = "Norwegian" OR movie_info.info = "German" OR movie_info.info = "USA" OR movie_info.info = "American") AND movie_info_idx.info > "6.0" AND title.production_year > 2010 AND (title.title LIKE "%murder%" OR title.title LIKE "%Murder%" OR title.title LIKE "%Mord%") AND kind_type.id = title.kind_id AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_info_idx.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id;

-- Query 17
-- 14c
SELECT movie_info_idx.info, title.title FROM info_type, info_type, keyword, kind_type, movie_info, movie_info_idx, movie_keyword, title WHERE info_type.info = "countries" AND info_type.info = "rating" AND (keyword.keyword = "murder" OR keyword.keyword = "murder-in-title" OR keyword.keyword = "blood" OR keyword.keyword = "violence") AND (kind_type.kind = "movie" OR kind_type.kind = "episode") AND (movie_info.info = "Sweden" OR movie_info.info = "Norway" OR movie_info.info = "Germany" OR movie_info.info = "Denmark" OR movie_info.info = "Swedish" OR movie_info.info = "Danish" OR movie_info.info = "Norwegian" OR movie_info.info = "German" OR movie_info.info = "USA" OR movie_info.info = "American") AND movie_info_idx.info < "8.5" AND title.production_year > 2005 AND kind_type.id = title.kind_id AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_info_idx.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id;

-- Query 18
-- 15a
SELECT movie_info.info, title.title FROM aka_title, company_name, company_type, info_type, keyword, movie_companies, movie_info, movie_keyword, title WHERE company_name.country_code = "[us]" AND info_type.info = "release dates" AND movie_companies.note LIKE "%(200%)%" AND movie_companies.note LIKE "%(worldwide)%" AND movie_info.note LIKE "%internet%" AND movie_info.info LIKE "USA:% 200%" AND title.production_year > 2000 AND title.id = aka_title.movie_id AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_companies.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = aka_title.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_info.movie_id = aka_title.movie_id AND movie_companies.movie_id = aka_title.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info.info_type_id AND company_name.id = movie_companies.company_id AND company_type.id = movie_companies.company_type_id;

-- Query 19
-- 15b
SELECT movie_info.info, title.title FROM aka_title, company_name, company_type, info_type, keyword, movie_companies, movie_info, movie_keyword, title WHERE company_name.country_code = "[us]" AND company_name.name = "YouTube" AND info_type.info = "release dates" AND movie_companies.note LIKE "%(200%)%" AND movie_companies.note LIKE "%(worldwide)%" AND movie_info.note LIKE "%internet%" AND movie_info.info LIKE "USA:% 200%" AND title.production_year >= 2005 AND title.production_year <= 2010 AND title.id = aka_title.movie_id AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_companies.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = aka_title.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_info.movie_id = aka_title.movie_id AND movie_companies.movie_id = aka_title.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info.info_type_id AND company_name.id = movie_companies.company_id AND company_type.id = movie_companies.company_type_id;

-- Query 20
-- 15c
SELECT movie_info.info, title.title FROM aka_title, company_name, company_type, info_type, keyword, movie_companies, movie_info, movie_keyword, title WHERE company_name.country_code = "[us]" AND info_type.info = "release dates" AND movie_info.note LIKE "%internet%" AND (movie_info.info LIKE "USA:% 199%" OR movie_info.info LIKE "USA:% 200%") AND title.production_year > 1990 AND title.id = aka_title.movie_id AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_companies.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = aka_title.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_info.movie_id = aka_title.movie_id AND movie_companies.movie_id = aka_title.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info.info_type_id AND company_name.id = movie_companies.company_id AND company_type.id = movie_companies.company_type_id;

-- Query 21
-- 15d
SELECT aka_title.title, title.title FROM aka_title, company_name, company_type, info_type, keyword, movie_companies, movie_info, movie_keyword, title WHERE company_name.country_code = "[us]" AND info_type.info = "release dates" AND movie_info.note LIKE "%internet%" AND title.production_year > 1990 AND title.id = aka_title.movie_id AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_companies.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = aka_title.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_info.movie_id = aka_title.movie_id AND movie_companies.movie_id = aka_title.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info.info_type_id AND company_name.id = movie_companies.company_id AND company_type.id = movie_companies.company_type_id;

-- Query 22
-- 16a
SELECT aka_name.name, title.title FROM aka_name, cast_info, company_name, keyword, movie_companies, movie_keyword, name, title WHERE company_name.country_code ="[us]" AND keyword.keyword ="character-name-in-title" AND title.episode_nr >= 50 AND title.episode_nr < 100 AND aka_name.person_id = name.id AND name.id = cast_info.person_id AND cast_info.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_id = company_name.id AND aka_name.person_id = cast_info.person_id AND cast_info.movie_id = movie_companies.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND movie_companies.movie_id = movie_keyword.movie_id;

-- Query 23
-- 16b
SELECT aka_name.name, title.title FROM aka_name, cast_info, company_name, keyword, movie_companies, movie_keyword, name, title WHERE company_name.country_code ="[us]" AND keyword.keyword ="character-name-in-title" AND aka_name.person_id = name.id AND name.id = cast_info.person_id AND cast_info.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_id = company_name.id AND aka_name.person_id = cast_info.person_id AND cast_info.movie_id = movie_companies.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND movie_companies.movie_id = movie_keyword.movie_id;

-- Query 24
-- 16c
SELECT aka_name.name, title.title FROM aka_name, cast_info, company_name, keyword, movie_companies, movie_keyword, name, title WHERE company_name.country_code ="[us]" AND keyword.keyword ="character-name-in-title" AND title.episode_nr < 100 AND aka_name.person_id = name.id AND name.id = cast_info.person_id AND cast_info.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_id = company_name.id AND aka_name.person_id = cast_info.person_id AND cast_info.movie_id = movie_companies.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND movie_companies.movie_id = movie_keyword.movie_id;

-- Query 25
-- 16d
SELECT aka_name.name, title.title FROM aka_name, cast_info, company_name, keyword, movie_companies, movie_keyword, name, title WHERE company_name.country_code ="[us]" AND keyword.keyword ="character-name-in-title" AND title.episode_nr >= 5 AND title.episode_nr < 100 AND aka_name.person_id = name.id AND name.id = cast_info.person_id AND cast_info.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_id = company_name.id AND aka_name.person_id = cast_info.person_id AND cast_info.movie_id = movie_companies.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND movie_companies.movie_id = movie_keyword.movie_id;

-- Query 26
-- 17a
SELECT name.name, name.name FROM cast_info, company_name, keyword, movie_companies, movie_keyword, name, title WHERE company_name.country_code ="[us]" AND keyword.keyword ="character-name-in-title" AND name.name LIKE "B%" AND name.id = cast_info.person_id AND cast_info.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_id = company_name.id AND cast_info.movie_id = movie_companies.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND movie_companies.movie_id = movie_keyword.movie_id;

-- Query 27
-- 17b
SELECT name.name, name.name FROM cast_info, company_name, keyword, movie_companies, movie_keyword, name, title WHERE keyword.keyword ="character-name-in-title" AND name.name LIKE "Z%" AND name.id = cast_info.person_id AND cast_info.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_id = company_name.id AND cast_info.movie_id = movie_companies.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND movie_companies.movie_id = movie_keyword.movie_id;

-- Query 28
-- 17c
SELECT name.name, name.name FROM cast_info, company_name, keyword, movie_companies, movie_keyword, name, title WHERE keyword.keyword ="character-name-in-title" AND name.name LIKE "X%" AND name.id = cast_info.person_id AND cast_info.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_id = company_name.id AND cast_info.movie_id = movie_companies.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND movie_companies.movie_id = movie_keyword.movie_id;

-- Query 29
-- 17d
SELECT name.name FROM cast_info, company_name, keyword, movie_companies, movie_keyword, name, title WHERE keyword.keyword ="character-name-in-title" AND name.name LIKE "%Bert%" AND name.id = cast_info.person_id AND cast_info.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_id = company_name.id AND cast_info.movie_id = movie_companies.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND movie_companies.movie_id = movie_keyword.movie_id;

-- Query 30
-- 17e
SELECT name.name FROM cast_info, company_name, keyword, movie_companies, movie_keyword, name, title WHERE company_name.country_code ="[us]" AND keyword.keyword ="character-name-in-title" AND name.id = cast_info.person_id AND cast_info.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_id = company_name.id AND cast_info.movie_id = movie_companies.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND movie_companies.movie_id = movie_keyword.movie_id;

-- Query 31
-- 17f
SELECT name.name FROM cast_info, company_name, keyword, movie_companies, movie_keyword, name, title WHERE keyword.keyword ="character-name-in-title" AND name.name LIKE "%B%" AND name.id = cast_info.person_id AND cast_info.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_id = company_name.id AND cast_info.movie_id = movie_companies.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND movie_companies.movie_id = movie_keyword.movie_id;

-- Query 32
-- 18a
SELECT movie_info.info, movie_info_idx.info, title.title FROM cast_info, info_type, info_type, movie_info, movie_info_idx, name, title WHERE (cast_info.note = "(producer)" OR cast_info.note = "(executive producer)") AND info_type.info = "budget" AND info_type.info = "votes" AND name.gender = "m" AND name.name LIKE "%Tim%" AND title.id = movie_info.movie_id AND title.id = movie_info_idx.movie_id AND title.id = cast_info.movie_id AND cast_info.movie_id = movie_info.movie_id AND cast_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND name.id = cast_info.person_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id;

-- Query 33
-- 18b
SELECT movie_info.info, movie_info_idx.info, title.title FROM cast_info, info_type, info_type, movie_info, movie_info_idx, name, title WHERE (cast_info.note = "(writer)" OR cast_info.note = "(head writer)" OR cast_info.note = "(written by)" OR cast_info.note = "(story)" OR cast_info.note = "(story editor)") AND info_type.info = "genres" AND info_type.info = "rating" AND (movie_info.info = "Horror" OR movie_info.info = "Thriller") AND movie_info_idx.info > "8.0" AND name.gender = "f" AND title.production_year >= 2008 AND title.production_year <= 2014 AND title.id = movie_info.movie_id AND title.id = movie_info_idx.movie_id AND title.id = cast_info.movie_id AND cast_info.movie_id = movie_info.movie_id AND cast_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND name.id = cast_info.person_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id;

-- Query 34
-- 18c
SELECT movie_info.info, movie_info_idx.info, title.title FROM cast_info, info_type, info_type, movie_info, movie_info_idx, name, title WHERE (cast_info.note = "(writer)" OR cast_info.note = "(head writer)" OR cast_info.note = "(written by)" OR cast_info.note = "(story)" OR cast_info.note = "(story editor)") AND info_type.info = "genres" AND info_type.info = "votes" AND (movie_info.info = "Horror" OR movie_info.info = "Action" OR movie_info.info = "Sci-Fi" OR movie_info.info = "Thriller" OR movie_info.info = "Crime" OR movie_info.info = "War") AND name.gender = "m" AND title.id = movie_info.movie_id AND title.id = movie_info_idx.movie_id AND title.id = cast_info.movie_id AND cast_info.movie_id = movie_info.movie_id AND cast_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND name.id = cast_info.person_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id;

-- Query 35
-- 19a
SELECT name.name, title.title FROM aka_name, char_name, cast_info, company_name, info_type, movie_companies, movie_info, name, role_type, title WHERE (cast_info.note = "(voice)" OR cast_info.note = "(voice: Japanese version)" OR cast_info.note = "(voice) (uncredited)" OR cast_info.note = "(voice: English version)") AND company_name.country_code ="[us]" AND info_type.info = "release dates" AND (movie_companies.note LIKE "%(USA)%" OR movie_companies.note LIKE "%(worldwide)%") AND (movie_info.info LIKE "Japan:%200%" OR movie_info.info LIKE "USA:%200%") AND name.gender ="f" AND name.name LIKE "%Ang%" AND role_type.role ="actress" AND title.production_year >= 2005 AND title.production_year <= 2009 AND title.id = movie_info.movie_id AND title.id = movie_companies.movie_id AND title.id = cast_info.movie_id AND movie_companies.movie_id = cast_info.movie_id AND movie_companies.movie_id = movie_info.movie_id AND movie_info.movie_id = cast_info.movie_id AND company_name.id = movie_companies.company_id AND info_type.id = movie_info.info_type_id AND name.id = cast_info.person_id AND role_type.id = cast_info.role_id AND name.id = aka_name.person_id AND cast_info.person_id = aka_name.person_id AND char_name.id = cast_info.person_role_id;

-- Query 36
-- 19b
SELECT name.name, title.title FROM aka_name, char_name, cast_info, company_name, info_type, movie_companies, movie_info, name, role_type, title WHERE cast_info.note = "(voice)" AND company_name.country_code ="[us]" AND info_type.info = "release dates" AND movie_companies.note LIKE "%(200%)%" AND (movie_companies.note LIKE "%(USA)%" OR movie_companies.note LIKE "%(worldwide)%") AND (movie_info.info LIKE "Japan:%2007%" OR movie_info.info LIKE "USA:%2008%") AND name.gender ="f" AND name.name LIKE "%Angel%" AND role_type.role ="actress" AND title.production_year >= 2007 AND title.production_year <= 2008 AND title.title LIKE "%Kung%Fu%Panda%" AND title.id = movie_info.movie_id AND title.id = movie_companies.movie_id AND title.id = cast_info.movie_id AND movie_companies.movie_id = cast_info.movie_id AND movie_companies.movie_id = movie_info.movie_id AND movie_info.movie_id = cast_info.movie_id AND company_name.id = movie_companies.company_id AND info_type.id = movie_info.info_type_id AND name.id = cast_info.person_id AND role_type.id = cast_info.role_id AND name.id = aka_name.person_id AND cast_info.person_id = aka_name.person_id AND char_name.id = cast_info.person_role_id;

-- Query 37
-- 19c
SELECT name.name, title.title FROM aka_name, char_name, cast_info, company_name, info_type, movie_companies, movie_info, name, role_type, title WHERE (cast_info.note = "(voice)" OR cast_info.note = "(voice: Japanese version)" OR cast_info.note = "(voice) (uncredited)" OR cast_info.note = "(voice: English version)") AND company_name.country_code ="[us]" AND info_type.info = "release dates" AND (movie_info.info LIKE "Japan:%200%" OR movie_info.info LIKE "USA:%200%") AND name.gender ="f" AND name.name LIKE "%An%" AND role_type.role ="actress" AND title.production_year > 2000 AND title.id = movie_info.movie_id AND title.id = movie_companies.movie_id AND title.id = cast_info.movie_id AND movie_companies.movie_id = cast_info.movie_id AND movie_companies.movie_id = movie_info.movie_id AND movie_info.movie_id = cast_info.movie_id AND company_name.id = movie_companies.company_id AND info_type.id = movie_info.info_type_id AND name.id = cast_info.person_id AND role_type.id = cast_info.role_id AND name.id = aka_name.person_id AND cast_info.person_id = aka_name.person_id AND char_name.id = cast_info.person_role_id;

-- Query 38
-- 19d
SELECT name.name, title.title FROM aka_name, char_name, cast_info, company_name, info_type, movie_companies, movie_info, name, role_type, title WHERE (cast_info.note = "(voice)" OR cast_info.note = "(voice: Japanese version)" OR cast_info.note = "(voice) (uncredited)" OR cast_info.note = "(voice: English version)") AND company_name.country_code ="[us]" AND info_type.info = "release dates" AND name.gender ="f" AND role_type.role ="actress" AND title.production_year > 2000 AND title.id = movie_info.movie_id AND title.id = movie_companies.movie_id AND title.id = cast_info.movie_id AND movie_companies.movie_id = cast_info.movie_id AND movie_companies.movie_id = movie_info.movie_id AND movie_info.movie_id = cast_info.movie_id AND company_name.id = movie_companies.company_id AND info_type.id = movie_info.info_type_id AND name.id = cast_info.person_id AND role_type.id = cast_info.role_id AND name.id = aka_name.person_id AND cast_info.person_id = aka_name.person_id AND char_name.id = cast_info.person_role_id;

-- Query 39
-- 1a
SELECT movie_companies.note, title.title, title.production_year FROM company_type, info_type, movie_companies, movie_info_idx, title WHERE company_type.kind = "production companies" AND info_type.info = "top 250 rank" AND (movie_companies.note LIKE "%(co-production)%" OR movie_companies.note LIKE "%(presents)%") AND company_type.id = movie_companies.company_type_id AND title.id = movie_companies.movie_id AND title.id = movie_info_idx.movie_id AND movie_companies.movie_id = movie_info_idx.movie_id AND info_type.id = movie_info_idx.info_type_id;

-- Query 40
-- 1b
SELECT movie_companies.note, title.title, title.production_year FROM company_type, info_type, movie_companies, movie_info_idx, title WHERE company_type.kind = "production companies" AND info_type.info = "bottom 10 rank" AND title.production_year >= 2005 AND title.production_year <= 2010 AND company_type.id = movie_companies.company_type_id AND title.id = movie_companies.movie_id AND title.id = movie_info_idx.movie_id AND movie_companies.movie_id = movie_info_idx.movie_id AND info_type.id = movie_info_idx.info_type_id;

-- Query 41
-- 1c
SELECT movie_companies.note, title.title, title.production_year FROM company_type, info_type, movie_companies, movie_info_idx, title WHERE company_type.kind = "production companies" AND info_type.info = "top 250 rank" AND (movie_companies.note LIKE "%(co-production)%") AND title.production_year >2010 AND company_type.id = movie_companies.company_type_id AND title.id = movie_companies.movie_id AND title.id = movie_info_idx.movie_id AND movie_companies.movie_id = movie_info_idx.movie_id AND info_type.id = movie_info_idx.info_type_id;

-- Query 42
-- 1d
SELECT movie_companies.note, title.title, title.production_year FROM company_type, info_type, movie_companies, movie_info_idx, title WHERE company_type.kind = "production companies" AND info_type.info = "bottom 10 rank" AND title.production_year >2000 AND company_type.id = movie_companies.company_type_id AND title.id = movie_companies.movie_id AND title.id = movie_info_idx.movie_id AND movie_companies.movie_id = movie_info_idx.movie_id AND info_type.id = movie_info_idx.info_type_id;

-- Query 43
-- 20a
SELECT title.title FROM complete_cast, comp_cast_type, comp_cast_type, char_name, cast_info, keyword, kind_type, movie_keyword, name, title WHERE comp_cast_type.kind = "cast" AND comp_cast_type.kind LIKE "%complete%" AND (char_name.name LIKE "%Tony%Stark%" OR char_name.name LIKE "%Iron%Man%") AND (keyword.keyword = "superhero" OR keyword.keyword = "sequel" OR keyword.keyword = "second-part" OR keyword.keyword = "marvel-comics" OR keyword.keyword = "based-on-comic" OR keyword.keyword = "tv-special" OR keyword.keyword = "fight" OR keyword.keyword = "violence") AND kind_type.kind = "movie" AND title.production_year > 1950 AND kind_type.id = title.kind_id AND title.id = movie_keyword.movie_id AND title.id = cast_info.movie_id AND title.id = complete_cast.movie_id AND movie_keyword.movie_id = cast_info.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND cast_info.movie_id = complete_cast.movie_id AND char_name.id = cast_info.person_role_id AND name.id = cast_info.person_id AND keyword.id = movie_keyword.keyword_id AND comp_cast_type.id = complete_cast.subject_id AND comp_cast_type.id = complete_cast.status_id;

-- Query 44
-- 20b
SELECT title.title FROM complete_cast, comp_cast_type, comp_cast_type, char_name, cast_info, keyword, kind_type, movie_keyword, name, title WHERE comp_cast_type.kind = "cast" AND comp_cast_type.kind LIKE "%complete%" AND (char_name.name LIKE "%Tony%Stark%" OR char_name.name LIKE "%Iron%Man%") AND (keyword.keyword = "superhero" OR keyword.keyword = "sequel" OR keyword.keyword = "second-part" OR keyword.keyword = "marvel-comics" OR keyword.keyword = "based-on-comic" OR keyword.keyword = "tv-special" OR keyword.keyword = "fight" OR keyword.keyword = "violence") AND kind_type.kind = "movie" AND name.name LIKE "%Downey%Robert%" AND title.production_year > 2000 AND kind_type.id = title.kind_id AND title.id = movie_keyword.movie_id AND title.id = cast_info.movie_id AND title.id = complete_cast.movie_id AND movie_keyword.movie_id = cast_info.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND cast_info.movie_id = complete_cast.movie_id AND char_name.id = cast_info.person_role_id AND name.id = cast_info.person_id AND keyword.id = movie_keyword.keyword_id AND comp_cast_type.id = complete_cast.subject_id AND comp_cast_type.id = complete_cast.status_id;

-- Query 45
-- 20c
SELECT name.name, title.title FROM complete_cast, comp_cast_type, comp_cast_type, char_name, cast_info, keyword, kind_type, movie_keyword, name, title WHERE comp_cast_type.kind = "cast" AND comp_cast_type.kind LIKE "%complete%" AND (char_name.name LIKE "%man%" OR char_name.name LIKE "%Man%") AND (keyword.keyword = "superhero" OR keyword.keyword = "marvel-comics" OR keyword.keyword = "based-on-comic" OR keyword.keyword = "tv-special" OR keyword.keyword = "fight" OR keyword.keyword = "violence" OR keyword.keyword = "magnet" OR keyword.keyword = "web" OR keyword.keyword = "claw" OR keyword.keyword = "laser") AND kind_type.kind = "movie" AND title.production_year > 2000 AND kind_type.id = title.kind_id AND title.id = movie_keyword.movie_id AND title.id = cast_info.movie_id AND title.id = complete_cast.movie_id AND movie_keyword.movie_id = cast_info.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND cast_info.movie_id = complete_cast.movie_id AND char_name.id = cast_info.person_role_id AND name.id = cast_info.person_id AND keyword.id = movie_keyword.keyword_id AND comp_cast_type.id = complete_cast.subject_id AND comp_cast_type.id = complete_cast.status_id;

-- Query 46
-- 21a
SELECT company_name.name, link_type.link, title.title FROM company_name, company_type, keyword, link_type, movie_companies, movie_info, movie_keyword, movie_link, title WHERE (company_name.name LIKE "%Film%" OR company_name.name LIKE "%Warner%") AND company_type.kind ="production companies" AND keyword.keyword ="sequel" AND link_type.link LIKE "%follow%" AND (movie_info.info = "Sweden" OR movie_info.info = "Norway" OR movie_info.info = "Germany" OR movie_info.info = "Denmark" OR movie_info.info = "Swedish" OR movie_info.info = "Denish" OR movie_info.info = "Norwegian" OR movie_info.info = "German") AND title.production_year >= 1950 AND title.production_year <= 2000 AND link_type.id = movie_link.link_type_id AND movie_link.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_type_id = company_type.id AND movie_companies.company_id = company_name.id AND movie_info.movie_id = title.id AND movie_link.movie_id = movie_keyword.movie_id AND movie_link.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_link.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_companies.movie_id = movie_info.movie_id;

-- Query 47
-- 21b
SELECT company_name.name, link_type.link, title.title FROM company_name, company_type, keyword, link_type, movie_companies, movie_info, movie_keyword, movie_link, title WHERE (company_name.name LIKE "%Film%" OR company_name.name LIKE "%Warner%") AND company_type.kind ="production companies" AND keyword.keyword ="sequel" AND link_type.link LIKE "%follow%" AND (movie_info.info = "Germany" OR movie_info.info = "German") AND title.production_year >= 2000 AND title.production_year <= 2010 AND link_type.id = movie_link.link_type_id AND movie_link.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_type_id = company_type.id AND movie_companies.company_id = company_name.id AND movie_info.movie_id = title.id AND movie_link.movie_id = movie_keyword.movie_id AND movie_link.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_link.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_companies.movie_id = movie_info.movie_id;

-- Query 48
-- 21c
SELECT company_name.name, link_type.link, title.title FROM company_name, company_type, keyword, link_type, movie_companies, movie_info, movie_keyword, movie_link, title WHERE (company_name.name LIKE "%Film%" OR company_name.name LIKE "%Warner%") AND company_type.kind ="production companies" AND keyword.keyword ="sequel" AND link_type.link LIKE "%follow%" AND (movie_info.info = "Sweden" OR movie_info.info = "Norway" OR movie_info.info = "Germany" OR movie_info.info = "Denmark" OR movie_info.info = "Swedish" OR movie_info.info = "Denish" OR movie_info.info = "Norwegian" OR movie_info.info = "German" OR movie_info.info = "English") AND title.production_year >= 1950 AND title.production_year <= 2010 AND link_type.id = movie_link.link_type_id AND movie_link.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_type_id = company_type.id AND movie_companies.company_id = company_name.id AND movie_info.movie_id = title.id AND movie_link.movie_id = movie_keyword.movie_id AND movie_link.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_link.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_companies.movie_id = movie_info.movie_id;

-- Query 49
-- 22a
SELECT company_name.name, movie_info_idx.info, title.title FROM company_name, company_type, info_type, info_type, keyword, kind_type, movie_companies, movie_info, movie_info_idx, movie_keyword, title WHERE info_type.info = "countries" AND info_type.info = "rating" AND (keyword.keyword = "murder" OR keyword.keyword = "murder-in-title" OR keyword.keyword = "blood" OR keyword.keyword = "violence") AND (kind_type.kind = "movie" OR kind_type.kind = "episode") AND movie_companies.note LIKE "%(200%)%" AND (movie_info.info = "Germany" OR movie_info.info = "German" OR movie_info.info = "USA" OR movie_info.info = "American") AND movie_info_idx.info < "7.0" AND title.production_year > 2008 AND kind_type.id = title.kind_id AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_info_idx.movie_id AND title.id = movie_companies.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_info_idx.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_companies.movie_id = movie_info_idx.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id AND company_type.id = movie_companies.company_type_id AND company_name.id = movie_companies.company_id;

-- Query 50
-- 22b
SELECT company_name.name, movie_info_idx.info, title.title FROM company_name, company_type, info_type, info_type, keyword, kind_type, movie_companies, movie_info, movie_info_idx, movie_keyword, title WHERE info_type.info = "countries" AND info_type.info = "rating" AND (keyword.keyword = "murder" OR keyword.keyword = "murder-in-title" OR keyword.keyword = "blood" OR keyword.keyword = "violence") AND (kind_type.kind = "movie" OR kind_type.kind = "episode") AND movie_companies.note LIKE "%(200%)%" AND (movie_info.info = "Germany" OR movie_info.info = "German" OR movie_info.info = "USA" OR movie_info.info = "American") AND movie_info_idx.info < "7.0" AND title.production_year > 2009 AND kind_type.id = title.kind_id AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_info_idx.movie_id AND title.id = movie_companies.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_info_idx.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_companies.movie_id = movie_info_idx.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id AND company_type.id = movie_companies.company_type_id AND company_name.id = movie_companies.company_id;

-- Query 51
-- 22c
SELECT company_name.name, movie_info_idx.info, title.title FROM company_name, company_type, info_type, info_type, keyword, kind_type, movie_companies, movie_info, movie_info_idx, movie_keyword, title WHERE info_type.info = "countries" AND info_type.info = "rating" AND (keyword.keyword = "murder" OR keyword.keyword = "murder-in-title" OR keyword.keyword = "blood" OR keyword.keyword = "violence") AND (kind_type.kind = "movie" OR kind_type.kind = "episode") AND movie_companies.note LIKE "%(200%)%" AND (movie_info.info = "Sweden" OR movie_info.info = "Norway" OR movie_info.info = "Germany" OR movie_info.info = "Denmark" OR movie_info.info = "Swedish" OR movie_info.info = "Danish" OR movie_info.info = "Norwegian" OR movie_info.info = "German" OR movie_info.info = "USA" OR movie_info.info = "American") AND movie_info_idx.info < "8.5" AND title.production_year > 2005 AND kind_type.id = title.kind_id AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_info_idx.movie_id AND title.id = movie_companies.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_info_idx.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_companies.movie_id = movie_info_idx.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id AND company_type.id = movie_companies.company_type_id AND company_name.id = movie_companies.company_id;

-- Query 52
-- 22d
SELECT company_name.name, movie_info_idx.info, title.title FROM company_name, company_type, info_type, info_type, keyword, kind_type, movie_companies, movie_info, movie_info_idx, movie_keyword, title WHERE info_type.info = "countries" AND info_type.info = "rating" AND (keyword.keyword = "murder" OR keyword.keyword = "murder-in-title" OR keyword.keyword = "blood" OR keyword.keyword = "violence") AND (kind_type.kind = "movie" OR kind_type.kind = "episode") AND (movie_info.info = "Sweden" OR movie_info.info = "Norway" OR movie_info.info = "Germany" OR movie_info.info = "Denmark" OR movie_info.info = "Swedish" OR movie_info.info = "Danish" OR movie_info.info = "Norwegian" OR movie_info.info = "German" OR movie_info.info = "USA" OR movie_info.info = "American") AND movie_info_idx.info < "8.5" AND title.production_year > 2005 AND kind_type.id = title.kind_id AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_info_idx.movie_id AND title.id = movie_companies.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_info_idx.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_companies.movie_id = movie_info_idx.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id AND company_type.id = movie_companies.company_type_id AND company_name.id = movie_companies.company_id;

-- Query 53
-- 23a
SELECT kind_type.kind, title.title FROM complete_cast, comp_cast_type, company_name, company_type, info_type, keyword, kind_type, movie_companies, movie_info, movie_keyword, title WHERE comp_cast_type.kind = "complete+verified" AND company_name.country_code = "[us]" AND info_type.info = "release dates" AND (kind_type.kind = "movie") AND movie_info.note LIKE "%internet%" AND (movie_info.info LIKE "USA:% 199%" OR movie_info.info LIKE "USA:% 200%") AND title.production_year > 2000 AND kind_type.id = title.kind_id AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_companies.movie_id AND title.id = complete_cast.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_info.movie_id = complete_cast.movie_id AND movie_companies.movie_id = complete_cast.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info.info_type_id AND company_name.id = movie_companies.company_id AND company_type.id = movie_companies.company_type_id AND comp_cast_type.id = complete_cast.status_id;

-- Query 54
-- 23b
SELECT kind_type.kind, title.title FROM complete_cast, comp_cast_type, company_name, company_type, info_type, keyword, kind_type, movie_companies, movie_info, movie_keyword, title WHERE comp_cast_type.kind = "complete+verified" AND company_name.country_code = "[us]" AND info_type.info = "release dates" AND (keyword.keyword = "nerd" OR keyword.keyword = "loner" OR keyword.keyword = "alienation" OR keyword.keyword = "dignity") AND (kind_type.kind = "movie") AND movie_info.note LIKE "%internet%" AND movie_info.info LIKE "USA:% 200%" AND title.production_year > 2000 AND kind_type.id = title.kind_id AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_companies.movie_id AND title.id = complete_cast.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_info.movie_id = complete_cast.movie_id AND movie_companies.movie_id = complete_cast.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info.info_type_id AND company_name.id = movie_companies.company_id AND company_type.id = movie_companies.company_type_id AND comp_cast_type.id = complete_cast.status_id;

-- Query 55
-- 23c
SELECT kind_type.kind, title.title FROM complete_cast, comp_cast_type, company_name, company_type, info_type, keyword, kind_type, movie_companies, movie_info, movie_keyword, title WHERE comp_cast_type.kind = "complete+verified" AND company_name.country_code = "[us]" AND info_type.info = "release dates" AND (kind_type.kind = "movie" OR kind_type.kind = "tv movie" OR kind_type.kind = "video movie" OR kind_type.kind = "video game") AND movie_info.note LIKE "%internet%" AND (movie_info.info LIKE "USA:% 199%" OR movie_info.info LIKE "USA:% 200%") AND title.production_year > 1990 AND kind_type.id = title.kind_id AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_companies.movie_id AND title.id = complete_cast.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_info.movie_id = complete_cast.movie_id AND movie_companies.movie_id = complete_cast.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info.info_type_id AND company_name.id = movie_companies.company_id AND company_type.id = movie_companies.company_type_id AND comp_cast_type.id = complete_cast.status_id;

-- Query 56
-- 24a
SELECT char_name.name, name.name, title.title FROM aka_name, char_name, cast_info, company_name, info_type, keyword, movie_companies, movie_info, movie_keyword, name, role_type, title WHERE (cast_info.note = "(voice)" OR cast_info.note = "(voice: Japanese version)" OR cast_info.note = "(voice) (uncredited)" OR cast_info.note = "(voice: English version)") AND company_name.country_code ="[us]" AND info_type.info = "release dates" AND (keyword.keyword = "hero" OR keyword.keyword = "martial-arts" OR keyword.keyword = "hand-to-hand-combat") AND (movie_info.info LIKE "Japan:%201%" OR movie_info.info LIKE "USA:%201%") AND name.gender ="f" AND name.name LIKE "%An%" AND role_type.role ="actress" AND title.production_year > 2010 AND title.id = movie_info.movie_id AND title.id = movie_companies.movie_id AND title.id = cast_info.movie_id AND title.id = movie_keyword.movie_id AND movie_companies.movie_id = cast_info.movie_id AND movie_companies.movie_id = movie_info.movie_id AND movie_companies.movie_id = movie_keyword.movie_id AND movie_info.movie_id = cast_info.movie_id AND movie_info.movie_id = movie_keyword.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND company_name.id = movie_companies.company_id AND info_type.id = movie_info.info_type_id AND name.id = cast_info.person_id AND role_type.id = cast_info.role_id AND name.id = aka_name.person_id AND cast_info.person_id = aka_name.person_id AND char_name.id = cast_info.person_role_id AND keyword.id = movie_keyword.keyword_id;

-- Query 57
-- 24b
SELECT char_name.name, name.name, title.title FROM aka_name, char_name, cast_info, company_name, info_type, keyword, movie_companies, movie_info, movie_keyword, name, role_type, title WHERE (cast_info.note = "(voice)" OR cast_info.note = "(voice: Japanese version)" OR cast_info.note = "(voice) (uncredited)" OR cast_info.note = "(voice: English version)") AND company_name.country_code ="[us]" AND company_name.name = "DreamWorks Animation" AND info_type.info = "release dates" AND (keyword.keyword = "hero" OR keyword.keyword = "martial-arts" OR keyword.keyword = "hand-to-hand-combat" OR keyword.keyword = "computer-animated-movie") AND (movie_info.info LIKE "Japan:%201%" OR movie_info.info LIKE "USA:%201%") AND name.gender ="f" AND name.name LIKE "%An%" AND role_type.role ="actress" AND title.production_year > 2010 AND title.title LIKE "Kung Fu Panda%" AND title.id = movie_info.movie_id AND title.id = movie_companies.movie_id AND title.id = cast_info.movie_id AND title.id = movie_keyword.movie_id AND movie_companies.movie_id = cast_info.movie_id AND movie_companies.movie_id = movie_info.movie_id AND movie_companies.movie_id = movie_keyword.movie_id AND movie_info.movie_id = cast_info.movie_id AND movie_info.movie_id = movie_keyword.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND company_name.id = movie_companies.company_id AND info_type.id = movie_info.info_type_id AND name.id = cast_info.person_id AND role_type.id = cast_info.role_id AND name.id = aka_name.person_id AND cast_info.person_id = aka_name.person_id AND char_name.id = cast_info.person_role_id AND keyword.id = movie_keyword.keyword_id;

-- Query 58
-- 25a
SELECT movie_info.info, movie_info_idx.info, name.name, title.title FROM cast_info, info_type, info_type, keyword, movie_info, movie_info_idx, movie_keyword, name, title WHERE (cast_info.note = "(writer)" OR cast_info.note = "(head writer)" OR cast_info.note = "(written by)" OR cast_info.note = "(story)" OR cast_info.note = "(story editor)") AND info_type.info = "genres" AND info_type.info = "votes" AND (keyword.keyword = "murder" OR keyword.keyword = "blood" OR keyword.keyword = "gore" OR keyword.keyword = "death" OR keyword.keyword = "female-nudity") AND movie_info.info = "Horror" AND name.gender = "m" AND title.id = movie_info.movie_id AND title.id = movie_info_idx.movie_id AND title.id = cast_info.movie_id AND title.id = movie_keyword.movie_id AND cast_info.movie_id = movie_info.movie_id AND cast_info.movie_id = movie_info_idx.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_keyword.movie_id AND movie_info_idx.movie_id = movie_keyword.movie_id AND name.id = cast_info.person_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id AND keyword.id = movie_keyword.keyword_id;

-- Query 59
-- 25b
SELECT movie_info.info, movie_info_idx.info, name.name, title.title FROM cast_info, info_type, info_type, keyword, movie_info, movie_info_idx, movie_keyword, name, title WHERE (cast_info.note = "(writer)" OR cast_info.note = "(head writer)" OR cast_info.note = "(written by)" OR cast_info.note = "(story)" OR cast_info.note = "(story editor)") AND info_type.info = "genres" AND info_type.info = "votes" AND (keyword.keyword = "murder" OR keyword.keyword = "blood" OR keyword.keyword = "gore" OR keyword.keyword = "death" OR keyword.keyword = "female-nudity") AND movie_info.info = "Horror" AND name.gender = "m" AND title.production_year > 2010 AND title.title LIKE "Vampire%" AND title.id = movie_info.movie_id AND title.id = movie_info_idx.movie_id AND title.id = cast_info.movie_id AND title.id = movie_keyword.movie_id AND cast_info.movie_id = movie_info.movie_id AND cast_info.movie_id = movie_info_idx.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_keyword.movie_id AND movie_info_idx.movie_id = movie_keyword.movie_id AND name.id = cast_info.person_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id AND keyword.id = movie_keyword.keyword_id;

-- Query 60
-- 25c
SELECT movie_info.info, movie_info_idx.info, name.name, title.title FROM cast_info, info_type, info_type, keyword, movie_info, movie_info_idx, movie_keyword, name, title WHERE (cast_info.note = "(writer)" OR cast_info.note = "(head writer)" OR cast_info.note = "(written by)" OR cast_info.note = "(story)" OR cast_info.note = "(story editor)") AND info_type.info = "genres" AND info_type.info = "votes" AND (keyword.keyword = "murder" OR keyword.keyword = "violence" OR keyword.keyword = "blood" OR keyword.keyword = "gore" OR keyword.keyword = "death" OR keyword.keyword = "female-nudity" OR keyword.keyword = "hospital") AND (movie_info.info = "Horror" OR movie_info.info = "Action" OR movie_info.info = "Sci-Fi" OR movie_info.info = "Thriller" OR movie_info.info = "Crime" OR movie_info.info = "War") AND name.gender = "m" AND title.id = movie_info.movie_id AND title.id = movie_info_idx.movie_id AND title.id = cast_info.movie_id AND title.id = movie_keyword.movie_id AND cast_info.movie_id = movie_info.movie_id AND cast_info.movie_id = movie_info_idx.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_keyword.movie_id AND movie_info_idx.movie_id = movie_keyword.movie_id AND name.id = cast_info.person_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id AND keyword.id = movie_keyword.keyword_id;

-- Query 61
-- 26a
SELECT char_name.name, movie_info_idx.info, name.name, title.title FROM complete_cast, comp_cast_type, comp_cast_type, char_name, cast_info, info_type, keyword, kind_type, movie_info_idx, movie_keyword, name, title WHERE comp_cast_type.kind = "cast" AND comp_cast_type.kind LIKE "%complete%" AND (char_name.name LIKE "%man%" OR char_name.name LIKE "%Man%") AND info_type.info = "rating" AND (keyword.keyword = "superhero" OR keyword.keyword = "marvel-comics" OR keyword.keyword = "based-on-comic" OR keyword.keyword = "tv-special" OR keyword.keyword = "fight" OR keyword.keyword = "violence" OR keyword.keyword = "magnet" OR keyword.keyword = "web" OR keyword.keyword = "claw" OR keyword.keyword = "laser") AND kind_type.kind = "movie" AND movie_info_idx.info > "7.0" AND title.production_year > 2000 AND kind_type.id = title.kind_id AND title.id = movie_keyword.movie_id AND title.id = cast_info.movie_id AND title.id = complete_cast.movie_id AND title.id = movie_info_idx.movie_id AND movie_keyword.movie_id = cast_info.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND movie_keyword.movie_id = movie_info_idx.movie_id AND cast_info.movie_id = complete_cast.movie_id AND cast_info.movie_id = movie_info_idx.movie_id AND complete_cast.movie_id = movie_info_idx.movie_id AND char_name.id = cast_info.person_role_id AND name.id = cast_info.person_id AND keyword.id = movie_keyword.keyword_id AND comp_cast_type.id = complete_cast.subject_id AND comp_cast_type.id = complete_cast.status_id AND info_type.id = movie_info_idx.info_type_id;

-- Query 62
-- 26b
SELECT char_name.name, movie_info_idx.info, title.title FROM complete_cast, comp_cast_type, comp_cast_type, char_name, cast_info, info_type, keyword, kind_type, movie_info_idx, movie_keyword, name, title WHERE comp_cast_type.kind = "cast" AND comp_cast_type.kind LIKE "%complete%" AND (char_name.name LIKE "%man%" OR char_name.name LIKE "%Man%") AND info_type.info = "rating" AND (keyword.keyword = "superhero" OR keyword.keyword = "marvel-comics" OR keyword.keyword = "based-on-comic" OR keyword.keyword = "fight") AND kind_type.kind = "movie" AND movie_info_idx.info > "8.0" AND title.production_year > 2005 AND kind_type.id = title.kind_id AND title.id = movie_keyword.movie_id AND title.id = cast_info.movie_id AND title.id = complete_cast.movie_id AND title.id = movie_info_idx.movie_id AND movie_keyword.movie_id = cast_info.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND movie_keyword.movie_id = movie_info_idx.movie_id AND cast_info.movie_id = complete_cast.movie_id AND cast_info.movie_id = movie_info_idx.movie_id AND complete_cast.movie_id = movie_info_idx.movie_id AND char_name.id = cast_info.person_role_id AND name.id = cast_info.person_id AND keyword.id = movie_keyword.keyword_id AND comp_cast_type.id = complete_cast.subject_id AND comp_cast_type.id = complete_cast.status_id AND info_type.id = movie_info_idx.info_type_id;

-- Query 63
-- 26c
SELECT char_name.name, movie_info_idx.info, title.title FROM complete_cast, comp_cast_type, comp_cast_type, char_name, cast_info, info_type, keyword, kind_type, movie_info_idx, movie_keyword, name, title WHERE comp_cast_type.kind = "cast" AND comp_cast_type.kind LIKE "%complete%" AND (char_name.name LIKE "%man%" OR char_name.name LIKE "%Man%") AND info_type.info = "rating" AND (keyword.keyword = "superhero" OR keyword.keyword = "marvel-comics" OR keyword.keyword = "based-on-comic" OR keyword.keyword = "tv-special" OR keyword.keyword = "fight" OR keyword.keyword = "violence" OR keyword.keyword = "magnet" OR keyword.keyword = "web" OR keyword.keyword = "claw" OR keyword.keyword = "laser") AND kind_type.kind = "movie" AND title.production_year > 2000 AND kind_type.id = title.kind_id AND title.id = movie_keyword.movie_id AND title.id = cast_info.movie_id AND title.id = complete_cast.movie_id AND title.id = movie_info_idx.movie_id AND movie_keyword.movie_id = cast_info.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND movie_keyword.movie_id = movie_info_idx.movie_id AND cast_info.movie_id = complete_cast.movie_id AND cast_info.movie_id = movie_info_idx.movie_id AND complete_cast.movie_id = movie_info_idx.movie_id AND char_name.id = cast_info.person_role_id AND name.id = cast_info.person_id AND keyword.id = movie_keyword.keyword_id AND comp_cast_type.id = complete_cast.subject_id AND comp_cast_type.id = complete_cast.status_id AND info_type.id = movie_info_idx.info_type_id;

-- Query 64
-- 27a
SELECT company_name.name, link_type.link, title.title FROM complete_cast, comp_cast_type, comp_cast_type, company_name, company_type, keyword, link_type, movie_companies, movie_info, movie_keyword, movie_link, title WHERE (comp_cast_type.kind = "cast" OR comp_cast_type.kind = "crew") AND comp_cast_type.kind = "complete" AND (company_name.name LIKE "%Film%" OR company_name.name LIKE "%Warner%") AND company_type.kind ="production companies" AND keyword.keyword ="sequel" AND link_type.link LIKE "%follow%" AND (movie_info.info = "Sweden" OR movie_info.info = "Germany" OR movie_info.info = "Swedish" OR movie_info.info = "German") AND title.production_year >= 1950 AND title.production_year <= 2000 AND link_type.id = movie_link.link_type_id AND movie_link.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_type_id = company_type.id AND movie_companies.company_id = company_name.id AND movie_info.movie_id = title.id AND title.id = complete_cast.movie_id AND comp_cast_type.id = complete_cast.subject_id AND comp_cast_type.id = complete_cast.status_id AND movie_link.movie_id = movie_keyword.movie_id AND movie_link.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_link.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_companies.movie_id = movie_info.movie_id AND movie_link.movie_id = complete_cast.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND movie_companies.movie_id = complete_cast.movie_id AND movie_info.movie_id = complete_cast.movie_id;

-- Query 65
-- 27b
SELECT company_name.name, link_type.link, title.title FROM complete_cast, comp_cast_type, comp_cast_type, company_name, company_type, keyword, link_type, movie_companies, movie_info, movie_keyword, movie_link, title WHERE (comp_cast_type.kind = "cast" OR comp_cast_type.kind = "crew") AND comp_cast_type.kind = "complete" AND (company_name.name LIKE "%Film%" OR company_name.name LIKE "%Warner%") AND company_type.kind ="production companies" AND keyword.keyword ="sequel" AND link_type.link LIKE "%follow%" AND (movie_info.info = "Sweden" OR movie_info.info = "Germany" OR movie_info.info = "Swedish" OR movie_info.info = "German") AND title.production_year = 1998 AND link_type.id = movie_link.link_type_id AND movie_link.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_type_id = company_type.id AND movie_companies.company_id = company_name.id AND movie_info.movie_id = title.id AND title.id = complete_cast.movie_id AND comp_cast_type.id = complete_cast.subject_id AND comp_cast_type.id = complete_cast.status_id AND movie_link.movie_id = movie_keyword.movie_id AND movie_link.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_link.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_companies.movie_id = movie_info.movie_id AND movie_link.movie_id = complete_cast.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND movie_companies.movie_id = complete_cast.movie_id AND movie_info.movie_id = complete_cast.movie_id;

-- Query 66
-- 27c
SELECT company_name.name, link_type.link, title.title FROM complete_cast, comp_cast_type, comp_cast_type, company_name, company_type, keyword, link_type, movie_companies, movie_info, movie_keyword, movie_link, title WHERE comp_cast_type.kind = "cast" AND comp_cast_type.kind LIKE "complete%" AND (company_name.name LIKE "%Film%" OR company_name.name LIKE "%Warner%") AND company_type.kind ="production companies" AND keyword.keyword ="sequel" AND link_type.link LIKE "%follow%" AND (movie_info.info = "Sweden" OR movie_info.info = "Norway" OR movie_info.info = "Germany" OR movie_info.info = "Denmark" OR movie_info.info = "Swedish" OR movie_info.info = "Denish" OR movie_info.info = "Norwegian" OR movie_info.info = "German" OR movie_info.info = "English") AND title.production_year >= 1950 AND title.production_year <= 2010 AND link_type.id = movie_link.link_type_id AND movie_link.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND title.id = movie_companies.movie_id AND movie_companies.company_type_id = company_type.id AND movie_companies.company_id = company_name.id AND movie_info.movie_id = title.id AND title.id = complete_cast.movie_id AND comp_cast_type.id = complete_cast.subject_id AND comp_cast_type.id = complete_cast.status_id AND movie_link.movie_id = movie_keyword.movie_id AND movie_link.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_link.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_companies.movie_id = movie_info.movie_id AND movie_link.movie_id = complete_cast.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND movie_companies.movie_id = complete_cast.movie_id AND movie_info.movie_id = complete_cast.movie_id;

-- Query 67
-- 28a
SELECT company_name.name, movie_info_idx.info, title.title FROM complete_cast, comp_cast_type, comp_cast_type, company_name, company_type, info_type, info_type, keyword, kind_type, movie_companies, movie_info, movie_info_idx, movie_keyword, title WHERE comp_cast_type.kind = "crew" AND info_type.info = "countries" AND info_type.info = "rating" AND (keyword.keyword = "murder" OR keyword.keyword = "murder-in-title" OR keyword.keyword = "blood" OR keyword.keyword = "violence") AND (kind_type.kind = "movie" OR kind_type.kind = "episode") AND movie_companies.note LIKE "%(200%)%" AND (movie_info.info = "Sweden" OR movie_info.info = "Norway" OR movie_info.info = "Germany" OR movie_info.info = "Denmark" OR movie_info.info = "Swedish" OR movie_info.info = "Danish" OR movie_info.info = "Norwegian" OR movie_info.info = "German" OR movie_info.info = "USA" OR movie_info.info = "American") AND movie_info_idx.info < "8.5" AND title.production_year > 2000 AND kind_type.id = title.kind_id AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_info_idx.movie_id AND title.id = movie_companies.movie_id AND title.id = complete_cast.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_info_idx.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_info.movie_id = complete_cast.movie_id AND movie_companies.movie_id = movie_info_idx.movie_id AND movie_companies.movie_id = complete_cast.movie_id AND movie_info_idx.movie_id = complete_cast.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id AND company_type.id = movie_companies.company_type_id AND company_name.id = movie_companies.company_id AND comp_cast_type.id = complete_cast.subject_id AND comp_cast_type.id = complete_cast.status_id;

-- Query 68
-- 28b
SELECT company_name.name, movie_info_idx.info, title.title FROM complete_cast, comp_cast_type, comp_cast_type, company_name, company_type, info_type, info_type, keyword, kind_type, movie_companies, movie_info, movie_info_idx, movie_keyword, title WHERE comp_cast_type.kind = "crew" AND info_type.info = "countries" AND info_type.info = "rating" AND (keyword.keyword = "murder" OR keyword.keyword = "murder-in-title" OR keyword.keyword = "blood" OR keyword.keyword = "violence") AND (kind_type.kind = "movie" OR kind_type.kind = "episode") AND movie_companies.note LIKE "%(200%)%" AND (movie_info.info = "Sweden" OR movie_info.info = "Germany" OR movie_info.info = "Swedish" OR movie_info.info = "German") AND movie_info_idx.info > "6.5" AND title.production_year > 2005 AND kind_type.id = title.kind_id AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_info_idx.movie_id AND title.id = movie_companies.movie_id AND title.id = complete_cast.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_info_idx.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_info.movie_id = complete_cast.movie_id AND movie_companies.movie_id = movie_info_idx.movie_id AND movie_companies.movie_id = complete_cast.movie_id AND movie_info_idx.movie_id = complete_cast.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id AND company_type.id = movie_companies.company_type_id AND company_name.id = movie_companies.company_id AND comp_cast_type.id = complete_cast.subject_id AND comp_cast_type.id = complete_cast.status_id;

-- Query 69
-- 28c
SELECT company_name.name, movie_info_idx.info, title.title FROM complete_cast, comp_cast_type, comp_cast_type, company_name, company_type, info_type, info_type, keyword, kind_type, movie_companies, movie_info, movie_info_idx, movie_keyword, title WHERE comp_cast_type.kind = "cast" AND comp_cast_type.kind = "complete" AND info_type.info = "countries" AND info_type.info = "rating" AND (keyword.keyword = "murder" OR keyword.keyword = "murder-in-title" OR keyword.keyword = "blood" OR keyword.keyword = "violence") AND (kind_type.kind = "movie" OR kind_type.kind = "episode") AND movie_companies.note LIKE "%(200%)%" AND (movie_info.info = "Sweden" OR movie_info.info = "Norway" OR movie_info.info = "Germany" OR movie_info.info = "Denmark" OR movie_info.info = "Swedish" OR movie_info.info = "Danish" OR movie_info.info = "Norwegian" OR movie_info.info = "German" OR movie_info.info = "USA" OR movie_info.info = "American") AND movie_info_idx.info < "8.5" AND title.production_year > 2005 AND kind_type.id = title.kind_id AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_info_idx.movie_id AND title.id = movie_companies.movie_id AND title.id = complete_cast.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND movie_keyword.movie_id = movie_info_idx.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_info.movie_id = complete_cast.movie_id AND movie_companies.movie_id = movie_info_idx.movie_id AND movie_companies.movie_id = complete_cast.movie_id AND movie_info_idx.movie_id = complete_cast.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id AND company_type.id = movie_companies.company_type_id AND company_name.id = movie_companies.company_id AND comp_cast_type.id = complete_cast.subject_id AND comp_cast_type.id = complete_cast.status_id;

-- Query 70
-- 29a
SELECT char_name.name, name.name, title.title FROM aka_name, complete_cast, comp_cast_type, comp_cast_type, char_name, cast_info, company_name, info_type, info_type, keyword, movie_companies, movie_info, movie_keyword, name, person_info, role_type, title WHERE comp_cast_type.kind ="cast" AND comp_cast_type.kind ="complete+verified" AND char_name.name = "Queen" AND (cast_info.note = "(voice)" OR cast_info.note = "(voice) (uncredited)" OR cast_info.note = "(voice: English version)") AND company_name.country_code ="[us]" AND info_type.info = "release dates" AND info_type.info = "trivia" AND keyword.keyword = "computer-animation" AND (movie_info.info LIKE "Japan:%200%" OR movie_info.info LIKE "USA:%200%") AND name.gender ="f" AND name.name LIKE "%An%" AND role_type.role ="actress" AND title.title = "Shrek 2" AND title.production_year >= 2000 AND title.production_year <= 2010 AND title.id = movie_info.movie_id AND title.id = movie_companies.movie_id AND title.id = cast_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = complete_cast.movie_id AND movie_companies.movie_id = cast_info.movie_id AND movie_companies.movie_id = movie_info.movie_id AND movie_companies.movie_id = movie_keyword.movie_id AND movie_companies.movie_id = complete_cast.movie_id AND movie_info.movie_id = cast_info.movie_id AND movie_info.movie_id = movie_keyword.movie_id AND movie_info.movie_id = complete_cast.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND cast_info.movie_id = complete_cast.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND company_name.id = movie_companies.company_id AND info_type.id = movie_info.info_type_id AND name.id = cast_info.person_id AND role_type.id = cast_info.role_id AND name.id = aka_name.person_id AND cast_info.person_id = aka_name.person_id AND char_name.id = cast_info.person_role_id AND name.id = person_info.person_id AND cast_info.person_id = person_info.person_id AND info_type.id = person_info.info_type_id AND keyword.id = movie_keyword.keyword_id AND comp_cast_type.id = complete_cast.subject_id AND comp_cast_type.id = complete_cast.status_id;

-- Query 71
-- 29b
SELECT char_name.name, name.name, title.title FROM aka_name, complete_cast, comp_cast_type, comp_cast_type, char_name, cast_info, company_name, info_type, info_type, keyword, movie_companies, movie_info, movie_keyword, name, person_info, role_type, title WHERE comp_cast_type.kind ="cast" AND comp_cast_type.kind ="complete+verified" AND char_name.name = "Queen" AND (cast_info.note = "(voice)" OR cast_info.note = "(voice) (uncredited)" OR cast_info.note = "(voice: English version)") AND company_name.country_code ="[us]" AND info_type.info = "release dates" AND info_type.info = "height" AND keyword.keyword = "computer-animation" AND movie_info.info LIKE "USA:%200%" AND name.gender ="f" AND name.name LIKE "%An%" AND role_type.role ="actress" AND title.title = "Shrek 2" AND title.production_year >= 2000 AND title.production_year <= 2005 AND title.id = movie_info.movie_id AND title.id = movie_companies.movie_id AND title.id = cast_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = complete_cast.movie_id AND movie_companies.movie_id = cast_info.movie_id AND movie_companies.movie_id = movie_info.movie_id AND movie_companies.movie_id = movie_keyword.movie_id AND movie_companies.movie_id = complete_cast.movie_id AND movie_info.movie_id = cast_info.movie_id AND movie_info.movie_id = movie_keyword.movie_id AND movie_info.movie_id = complete_cast.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND cast_info.movie_id = complete_cast.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND company_name.id = movie_companies.company_id AND info_type.id = movie_info.info_type_id AND name.id = cast_info.person_id AND role_type.id = cast_info.role_id AND name.id = aka_name.person_id AND cast_info.person_id = aka_name.person_id AND char_name.id = cast_info.person_role_id AND name.id = person_info.person_id AND cast_info.person_id = person_info.person_id AND info_type.id = person_info.info_type_id AND keyword.id = movie_keyword.keyword_id AND comp_cast_type.id = complete_cast.subject_id AND comp_cast_type.id = complete_cast.status_id;

-- Query 72
-- 29c
SELECT char_name.name, name.name, title.title FROM aka_name, complete_cast, comp_cast_type, comp_cast_type, char_name, cast_info, company_name, info_type, info_type, keyword, movie_companies, movie_info, movie_keyword, name, person_info, role_type, title WHERE comp_cast_type.kind ="cast" AND comp_cast_type.kind ="complete+verified" AND (cast_info.note = "(voice)" OR cast_info.note = "(voice: Japanese version)" OR cast_info.note = "(voice) (uncredited)" OR cast_info.note = "(voice: English version)") AND company_name.country_code ="[us]" AND info_type.info = "release dates" AND info_type.info = "trivia" AND keyword.keyword = "computer-animation" AND (movie_info.info LIKE "Japan:%200%" OR movie_info.info LIKE "USA:%200%") AND name.gender ="f" AND name.name LIKE "%An%" AND role_type.role ="actress" AND title.production_year >= 2000 AND title.production_year <= 2010 AND title.id = movie_info.movie_id AND title.id = movie_companies.movie_id AND title.id = cast_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = complete_cast.movie_id AND movie_companies.movie_id = cast_info.movie_id AND movie_companies.movie_id = movie_info.movie_id AND movie_companies.movie_id = movie_keyword.movie_id AND movie_companies.movie_id = complete_cast.movie_id AND movie_info.movie_id = cast_info.movie_id AND movie_info.movie_id = movie_keyword.movie_id AND movie_info.movie_id = complete_cast.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND cast_info.movie_id = complete_cast.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND company_name.id = movie_companies.company_id AND info_type.id = movie_info.info_type_id AND name.id = cast_info.person_id AND role_type.id = cast_info.role_id AND name.id = aka_name.person_id AND cast_info.person_id = aka_name.person_id AND char_name.id = cast_info.person_role_id AND name.id = person_info.person_id AND cast_info.person_id = person_info.person_id AND info_type.id = person_info.info_type_id AND keyword.id = movie_keyword.keyword_id AND comp_cast_type.id = complete_cast.subject_id AND comp_cast_type.id = complete_cast.status_id;

-- Query 73
-- 2a
SELECT title.title FROM company_name, keyword, movie_companies, movie_keyword, title WHERE company_name.country_code ="[de]" AND keyword.keyword ="character-name-in-title" AND company_name.id = movie_companies.company_id AND movie_companies.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND movie_companies.movie_id = movie_keyword.movie_id;

-- Query 74
-- 2b
SELECT title.title FROM company_name, keyword, movie_companies, movie_keyword, title WHERE company_name.country_code ="[nl]" AND keyword.keyword ="character-name-in-title" AND company_name.id = movie_companies.company_id AND movie_companies.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND movie_companies.movie_id = movie_keyword.movie_id;

-- Query 75
-- 2c
SELECT title.title FROM company_name, keyword, movie_companies, movie_keyword, title WHERE company_name.country_code ="[sm]" AND keyword.keyword ="character-name-in-title" AND company_name.id = movie_companies.company_id AND movie_companies.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND movie_companies.movie_id = movie_keyword.movie_id;

-- Query 76
-- 2d
SELECT title.title FROM company_name, keyword, movie_companies, movie_keyword, title WHERE company_name.country_code ="[us]" AND keyword.keyword ="character-name-in-title" AND company_name.id = movie_companies.company_id AND movie_companies.movie_id = title.id AND title.id = movie_keyword.movie_id AND movie_keyword.keyword_id = keyword.id AND movie_companies.movie_id = movie_keyword.movie_id;

-- Query 77
-- 30a
SELECT movie_info.info, movie_info_idx.info, name.name, title.title FROM complete_cast, comp_cast_type, comp_cast_type, cast_info, info_type, info_type, keyword, movie_info, movie_info_idx, movie_keyword, name, title WHERE (comp_cast_type.kind = "cast" OR comp_cast_type.kind = "crew") AND comp_cast_type.kind ="complete+verified" AND (cast_info.note = "(writer)" OR cast_info.note = "(head writer)" OR cast_info.note = "(written by)" OR cast_info.note = "(story)" OR cast_info.note = "(story editor)") AND info_type.info = "genres" AND info_type.info = "votes" AND (keyword.keyword = "murder" OR keyword.keyword = "violence" OR keyword.keyword = "blood" OR keyword.keyword = "gore" OR keyword.keyword = "death" OR keyword.keyword = "female-nudity" OR keyword.keyword = "hospital") AND (movie_info.info = "Horror" OR movie_info.info = "Thriller") AND name.gender = "m" AND title.production_year > 2000 AND title.id = movie_info.movie_id AND title.id = movie_info_idx.movie_id AND title.id = cast_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = complete_cast.movie_id AND cast_info.movie_id = movie_info.movie_id AND cast_info.movie_id = movie_info_idx.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND cast_info.movie_id = complete_cast.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_keyword.movie_id AND movie_info.movie_id = complete_cast.movie_id AND movie_info_idx.movie_id = movie_keyword.movie_id AND movie_info_idx.movie_id = complete_cast.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND name.id = cast_info.person_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id AND keyword.id = movie_keyword.keyword_id AND comp_cast_type.id = complete_cast.subject_id AND comp_cast_type.id = complete_cast.status_id;

-- Query 78
-- 30b
SELECT movie_info.info, movie_info_idx.info, name.name, title.title FROM complete_cast, comp_cast_type, comp_cast_type, cast_info, info_type, info_type, keyword, movie_info, movie_info_idx, movie_keyword, name, title WHERE (comp_cast_type.kind = "cast" OR comp_cast_type.kind = "crew") AND comp_cast_type.kind ="complete+verified" AND (cast_info.note = "(writer)" OR cast_info.note = "(head writer)" OR cast_info.note = "(written by)" OR cast_info.note = "(story)" OR cast_info.note = "(story editor)") AND info_type.info = "genres" AND info_type.info = "votes" AND (keyword.keyword = "murder" OR keyword.keyword = "violence" OR keyword.keyword = "blood" OR keyword.keyword = "gore" OR keyword.keyword = "death" OR keyword.keyword = "female-nudity" OR keyword.keyword = "hospital") AND (movie_info.info = "Horror" OR movie_info.info = "Thriller") AND name.gender = "m" AND title.production_year > 2000 AND (title.title LIKE "%Freddy%" OR title.title LIKE "%Jason%" OR title.title LIKE "Saw%") AND title.id = movie_info.movie_id AND title.id = movie_info_idx.movie_id AND title.id = cast_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = complete_cast.movie_id AND cast_info.movie_id = movie_info.movie_id AND cast_info.movie_id = movie_info_idx.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND cast_info.movie_id = complete_cast.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_keyword.movie_id AND movie_info.movie_id = complete_cast.movie_id AND movie_info_idx.movie_id = movie_keyword.movie_id AND movie_info_idx.movie_id = complete_cast.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND name.id = cast_info.person_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id AND keyword.id = movie_keyword.keyword_id AND comp_cast_type.id = complete_cast.subject_id AND comp_cast_type.id = complete_cast.status_id;

-- Query 79
-- 30c
SELECT movie_info.info, movie_info_idx.info, name.name, title.title FROM complete_cast, comp_cast_type, comp_cast_type, cast_info, info_type, info_type, keyword, movie_info, movie_info_idx, movie_keyword, name, title WHERE comp_cast_type.kind = "cast" AND comp_cast_type.kind ="complete+verified" AND (cast_info.note = "(writer)" OR cast_info.note = "(head writer)" OR cast_info.note = "(written by)" OR cast_info.note = "(story)" OR cast_info.note = "(story editor)") AND info_type.info = "genres" AND info_type.info = "votes" AND (keyword.keyword = "murder" OR keyword.keyword = "violence" OR keyword.keyword = "blood" OR keyword.keyword = "gore" OR keyword.keyword = "death" OR keyword.keyword = "female-nudity" OR keyword.keyword = "hospital") AND (movie_info.info = "Horror" OR movie_info.info = "Action" OR movie_info.info = "Sci-Fi" OR movie_info.info = "Thriller" OR movie_info.info = "Crime" OR movie_info.info = "War") AND name.gender = "m" AND title.id = movie_info.movie_id AND title.id = movie_info_idx.movie_id AND title.id = cast_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = complete_cast.movie_id AND cast_info.movie_id = movie_info.movie_id AND cast_info.movie_id = movie_info_idx.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND cast_info.movie_id = complete_cast.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_keyword.movie_id AND movie_info.movie_id = complete_cast.movie_id AND movie_info_idx.movie_id = movie_keyword.movie_id AND movie_info_idx.movie_id = complete_cast.movie_id AND movie_keyword.movie_id = complete_cast.movie_id AND name.id = cast_info.person_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id AND keyword.id = movie_keyword.keyword_id AND comp_cast_type.id = complete_cast.subject_id AND comp_cast_type.id = complete_cast.status_id;

-- Query 80
-- 31a
SELECT movie_info.info, movie_info_idx.info, name.name, title.title FROM cast_info, company_name, info_type, info_type, keyword, movie_companies, movie_info, movie_info_idx, movie_keyword, name, title WHERE (cast_info.note = "(writer)" OR cast_info.note = "(head writer)" OR cast_info.note = "(written by)" OR cast_info.note = "(story)" OR cast_info.note = "(story editor)") AND company_name.name LIKE "Lionsgate%" AND info_type.info = "genres" AND info_type.info = "votes" AND (keyword.keyword = "murder" OR keyword.keyword = "violence" OR keyword.keyword = "blood" OR keyword.keyword = "gore" OR keyword.keyword = "death" OR keyword.keyword = "female-nudity" OR keyword.keyword = "hospital") AND (movie_info.info = "Horror" OR movie_info.info = "Thriller") AND name.gender = "m" AND title.id = movie_info.movie_id AND title.id = movie_info_idx.movie_id AND title.id = cast_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_companies.movie_id AND cast_info.movie_id = movie_info.movie_id AND cast_info.movie_id = movie_info_idx.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND cast_info.movie_id = movie_companies.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_keyword.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_info_idx.movie_id = movie_keyword.movie_id AND movie_info_idx.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND name.id = cast_info.person_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id AND keyword.id = movie_keyword.keyword_id AND company_name.id = movie_companies.company_id;

-- Query 81
-- 31b
SELECT movie_info.info, movie_info_idx.info, name.name, title.title FROM cast_info, company_name, info_type, info_type, keyword, movie_companies, movie_info, movie_info_idx, movie_keyword, name, title WHERE (cast_info.note = "(writer)" OR cast_info.note = "(head writer)" OR cast_info.note = "(written by)" OR cast_info.note = "(story)" OR cast_info.note = "(story editor)") AND company_name.name LIKE "Lionsgate%" AND info_type.info = "genres" AND info_type.info = "votes" AND (keyword.keyword = "murder" OR keyword.keyword = "violence" OR keyword.keyword = "blood" OR keyword.keyword = "gore" OR keyword.keyword = "death" OR keyword.keyword = "female-nudity" OR keyword.keyword = "hospital") AND movie_companies.note LIKE "%(Blu-ray)%" AND (movie_info.info = "Horror" OR movie_info.info = "Thriller") AND name.gender = "m" AND title.production_year > 2000 AND (title.title LIKE "%Freddy%" OR title.title LIKE "%Jason%" OR title.title LIKE "Saw%") AND title.id = movie_info.movie_id AND title.id = movie_info_idx.movie_id AND title.id = cast_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_companies.movie_id AND cast_info.movie_id = movie_info.movie_id AND cast_info.movie_id = movie_info_idx.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND cast_info.movie_id = movie_companies.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_keyword.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_info_idx.movie_id = movie_keyword.movie_id AND movie_info_idx.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND name.id = cast_info.person_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id AND keyword.id = movie_keyword.keyword_id AND company_name.id = movie_companies.company_id;

-- Query 82
-- 31c
SELECT movie_info.info, movie_info_idx.info, name.name, title.title FROM cast_info, company_name, info_type, info_type, keyword, movie_companies, movie_info, movie_info_idx, movie_keyword, name, title WHERE (cast_info.note = "(writer)" OR cast_info.note = "(head writer)" OR cast_info.note = "(written by)" OR cast_info.note = "(story)" OR cast_info.note = "(story editor)") AND company_name.name LIKE "Lionsgate%" AND info_type.info = "genres" AND info_type.info = "votes" AND (keyword.keyword = "murder" OR keyword.keyword = "violence" OR keyword.keyword = "blood" OR keyword.keyword = "gore" OR keyword.keyword = "death" OR keyword.keyword = "female-nudity" OR keyword.keyword = "hospital") AND (movie_info.info = "Horror" OR movie_info.info = "Action" OR movie_info.info = "Sci-Fi" OR movie_info.info = "Thriller" OR movie_info.info = "Crime" OR movie_info.info = "War") AND title.id = movie_info.movie_id AND title.id = movie_info_idx.movie_id AND title.id = cast_info.movie_id AND title.id = movie_keyword.movie_id AND title.id = movie_companies.movie_id AND cast_info.movie_id = movie_info.movie_id AND cast_info.movie_id = movie_info_idx.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND cast_info.movie_id = movie_companies.movie_id AND movie_info.movie_id = movie_info_idx.movie_id AND movie_info.movie_id = movie_keyword.movie_id AND movie_info.movie_id = movie_companies.movie_id AND movie_info_idx.movie_id = movie_keyword.movie_id AND movie_info_idx.movie_id = movie_companies.movie_id AND movie_keyword.movie_id = movie_companies.movie_id AND name.id = cast_info.person_id AND info_type.id = movie_info.info_type_id AND info_type.id = movie_info_idx.info_type_id AND keyword.id = movie_keyword.keyword_id AND company_name.id = movie_companies.company_id;

-- Query 83
-- 32a
SELECT link_type.link, title.title, title.title FROM keyword, link_type, movie_keyword, movie_link, title, title WHERE keyword.keyword ="10,000-mile-club" AND movie_keyword.keyword_id = keyword.id AND title.id = movie_keyword.movie_id AND movie_link.movie_id = title.id AND movie_link.linked_movie_id = title.id AND link_type.id = movie_link.link_type_id AND movie_keyword.movie_id = title.id;

-- Query 84
-- 32b
SELECT link_type.link, title.title, title.title FROM keyword, link_type, movie_keyword, movie_link, title, title WHERE keyword.keyword ="character-name-in-title" AND movie_keyword.keyword_id = keyword.id AND title.id = movie_keyword.movie_id AND movie_link.movie_id = title.id AND movie_link.linked_movie_id = title.id AND link_type.id = movie_link.link_type_id AND movie_keyword.movie_id = title.id;

-- Query 85
-- 33a
SELECT company_name.name, company_name.name, movie_info_idx.info, movie_info_idx.info, title.title, title.title FROM company_name, company_name, info_type, info_type, kind_type, kind_type, link_type, movie_companies, movie_companies, movie_info_idx, movie_info_idx, movie_link, title, title WHERE company_name.country_code = "[us]" AND info_type.info = "rating" AND info_type.info = "rating" AND (kind_type.kind = "tv series") AND (kind_type.kind = "tv series") AND (link_type.link = "sequel" OR link_type.link = "follows" OR link_type.link = "followed by") AND movie_info_idx.info < "3.0" AND title.production_year >= 2005 AND title.production_year <= 2008 AND link_type.id = movie_link.link_type_id AND title.id = movie_link.movie_id AND title.id = movie_link.linked_movie_id AND info_type.id = movie_info_idx.info_type_id AND title.id = movie_info_idx.movie_id AND kind_type.id = title.kind_id AND company_name.id = movie_companies.company_id AND title.id = movie_companies.movie_id AND movie_link.movie_id = movie_info_idx.movie_id AND movie_link.movie_id = movie_companies.movie_id AND movie_info_idx.movie_id = movie_companies.movie_id AND info_type.id = movie_info_idx.info_type_id AND title.id = movie_info_idx.movie_id AND kind_type.id = title.kind_id AND company_name.id = movie_companies.company_id AND title.id = movie_companies.movie_id AND movie_link.linked_movie_id = movie_info_idx.movie_id AND movie_link.linked_movie_id = movie_companies.movie_id AND movie_info_idx.movie_id = movie_companies.movie_id;

-- Query 86
-- 33b
SELECT company_name.name, company_name.name, movie_info_idx.info, movie_info_idx.info, title.title, title.title FROM company_name, company_name, info_type, info_type, kind_type, kind_type, link_type, movie_companies, movie_companies, movie_info_idx, movie_info_idx, movie_link, title, title WHERE company_name.country_code = "[nl]" AND info_type.info = "rating" AND info_type.info = "rating" AND (kind_type.kind = "tv series") AND (kind_type.kind = "tv series") AND link_type.link LIKE "%follow%" AND movie_info_idx.info < "3.0" AND title.production_year = 2007 AND link_type.id = movie_link.link_type_id AND title.id = movie_link.movie_id AND title.id = movie_link.linked_movie_id AND info_type.id = movie_info_idx.info_type_id AND title.id = movie_info_idx.movie_id AND kind_type.id = title.kind_id AND company_name.id = movie_companies.company_id AND title.id = movie_companies.movie_id AND movie_link.movie_id = movie_info_idx.movie_id AND movie_link.movie_id = movie_companies.movie_id AND movie_info_idx.movie_id = movie_companies.movie_id AND info_type.id = movie_info_idx.info_type_id AND title.id = movie_info_idx.movie_id AND kind_type.id = title.kind_id AND company_name.id = movie_companies.company_id AND title.id = movie_companies.movie_id AND movie_link.linked_movie_id = movie_info_idx.movie_id AND movie_link.linked_movie_id = movie_companies.movie_id AND movie_info_idx.movie_id = movie_companies.movie_id;

-- Query 87
-- 33c
SELECT company_name.name, company_name.name, movie_info_idx.info, movie_info_idx.info, title.title, title.title FROM company_name, company_name, info_type, info_type, kind_type, kind_type, link_type, movie_companies, movie_companies, movie_info_idx, movie_info_idx, movie_link, title, title WHERE info_type.info = "rating" AND info_type.info = "rating" AND (kind_type.kind = "tv series" OR kind_type.kind = "episode") AND (kind_type.kind = "tv series" OR kind_type.kind = "episode") AND (link_type.link = "sequel" OR link_type.link = "follows" OR link_type.link = "followed by") AND movie_info_idx.info < "3.5" AND title.production_year >= 2000 AND title.production_year <= 2010 AND link_type.id = movie_link.link_type_id AND title.id = movie_link.movie_id AND title.id = movie_link.linked_movie_id AND info_type.id = movie_info_idx.info_type_id AND title.id = movie_info_idx.movie_id AND kind_type.id = title.kind_id AND company_name.id = movie_companies.company_id AND title.id = movie_companies.movie_id AND movie_link.movie_id = movie_info_idx.movie_id AND movie_link.movie_id = movie_companies.movie_id AND movie_info_idx.movie_id = movie_companies.movie_id AND info_type.id = movie_info_idx.info_type_id AND title.id = movie_info_idx.movie_id AND kind_type.id = title.kind_id AND company_name.id = movie_companies.company_id AND title.id = movie_companies.movie_id AND movie_link.linked_movie_id = movie_info_idx.movie_id AND movie_link.linked_movie_id = movie_companies.movie_id AND movie_info_idx.movie_id = movie_companies.movie_id;

-- Query 88
-- 3a
SELECT title.title FROM keyword, movie_info, movie_keyword, title WHERE keyword.keyword LIKE "%sequel%" AND (movie_info.info = "Sweden" OR movie_info.info = "Norway" OR movie_info.info = "Germany" OR movie_info.info = "Denmark" OR movie_info.info = "Swedish" OR movie_info.info = "Denish" OR movie_info.info = "Norwegian" OR movie_info.info = "German") AND title.production_year > 2005 AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND keyword.id = movie_keyword.keyword_id;

-- Query 89
-- 3b
SELECT title.title FROM keyword, movie_info, movie_keyword, title WHERE keyword.keyword LIKE "%sequel%" AND (movie_info.info = "Bulgaria") AND title.production_year > 2010 AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND keyword.id = movie_keyword.keyword_id;

-- Query 90
-- 3c
SELECT title.title FROM keyword, movie_info, movie_keyword, title WHERE keyword.keyword LIKE "%sequel%" AND (movie_info.info = "Sweden" OR movie_info.info = "Norway" OR movie_info.info = "Germany" OR movie_info.info = "Denmark" OR movie_info.info = "Swedish" OR movie_info.info = "Denish" OR movie_info.info = "Norwegian" OR movie_info.info = "German" OR movie_info.info = "USA" OR movie_info.info = "American") AND title.production_year > 1990 AND title.id = movie_info.movie_id AND title.id = movie_keyword.movie_id AND movie_keyword.movie_id = movie_info.movie_id AND keyword.id = movie_keyword.keyword_id;

-- Query 91
-- 4a
SELECT movie_info_idx.info, title.title FROM info_type, keyword, movie_info_idx, movie_keyword, title WHERE info_type.info ="rating" AND keyword.keyword LIKE "%sequel%" AND movie_info_idx.info > "5.0" AND title.production_year > 2005 AND title.id = movie_info_idx.movie_id AND title.id = movie_keyword.movie_id AND movie_keyword.movie_id = movie_info_idx.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info_idx.info_type_id;

-- Query 92
-- 4b
SELECT movie_info_idx.info, title.title FROM info_type, keyword, movie_info_idx, movie_keyword, title WHERE info_type.info ="rating" AND keyword.keyword LIKE "%sequel%" AND movie_info_idx.info > "9.0" AND title.production_year > 2010 AND title.id = movie_info_idx.movie_id AND title.id = movie_keyword.movie_id AND movie_keyword.movie_id = movie_info_idx.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info_idx.info_type_id;

-- Query 93
-- 4c
SELECT movie_info_idx.info, title.title FROM info_type, keyword, movie_info_idx, movie_keyword, title WHERE info_type.info ="rating" AND keyword.keyword LIKE "%sequel%" AND movie_info_idx.info > "2.0" AND title.production_year > 1990 AND title.id = movie_info_idx.movie_id AND title.id = movie_keyword.movie_id AND movie_keyword.movie_id = movie_info_idx.movie_id AND keyword.id = movie_keyword.keyword_id AND info_type.id = movie_info_idx.info_type_id;

-- Query 94
-- 5a
SELECT title.title FROM company_type, info_type, movie_companies, movie_info, title WHERE company_type.kind = "production companies" AND movie_companies.note LIKE "%(theatrical)%" AND movie_companies.note LIKE "%(France)%" AND (movie_info.info = "Sweden" OR movie_info.info = "Norway" OR movie_info.info = "Germany" OR movie_info.info = "Denmark" OR movie_info.info = "Swedish" OR movie_info.info = "Denish" OR movie_info.info = "Norwegian" OR movie_info.info = "German") AND title.production_year > 2005 AND title.id = movie_info.movie_id AND title.id = movie_companies.movie_id AND movie_companies.movie_id = movie_info.movie_id AND company_type.id = movie_companies.company_type_id AND info_type.id = movie_info.info_type_id;

-- Query 95
-- 5b
SELECT title.title FROM company_type, info_type, movie_companies, movie_info, title WHERE company_type.kind = "production companies" AND movie_companies.note LIKE "%(VHS)%" AND movie_companies.note LIKE "%(USA)%" AND movie_companies.note LIKE "%(1994)%" AND (movie_info.info = "USA" OR movie_info.info = "America") AND title.production_year > 2010 AND title.id = movie_info.movie_id AND title.id = movie_companies.movie_id AND movie_companies.movie_id = movie_info.movie_id AND company_type.id = movie_companies.company_type_id AND info_type.id = movie_info.info_type_id;

-- Query 96
-- 5c
SELECT title.title FROM company_type, info_type, movie_companies, movie_info, title WHERE company_type.kind = "production companies" AND movie_companies.note LIKE "%(USA)%" AND (movie_info.info = "Sweden" OR movie_info.info = "Norway" OR movie_info.info = "Germany" OR movie_info.info = "Denmark" OR movie_info.info = "Swedish" OR movie_info.info = "Denish" OR movie_info.info = "Norwegian" OR movie_info.info = "German" OR movie_info.info = "USA" OR movie_info.info = "American") AND title.production_year > 1990 AND title.id = movie_info.movie_id AND title.id = movie_companies.movie_id AND movie_companies.movie_id = movie_info.movie_id AND company_type.id = movie_companies.company_type_id AND info_type.id = movie_info.info_type_id;

-- Query 97
-- 6a
SELECT keyword.keyword, name.name, title.title FROM cast_info, keyword, movie_keyword, name, title WHERE keyword.keyword = "marvel-cinematic-universe" AND name.name LIKE "%Downey%Robert%" AND title.production_year > 2010 AND keyword.id = movie_keyword.keyword_id AND title.id = movie_keyword.movie_id AND title.id = cast_info.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND name.id = cast_info.person_id;

-- Query 98
-- 6b
SELECT keyword.keyword, name.name, title.title FROM cast_info, keyword, movie_keyword, name, title WHERE (keyword.keyword = "superhero" OR keyword.keyword = "sequel" OR keyword.keyword = "second-part" OR keyword.keyword = "marvel-comics" OR keyword.keyword = "based-on-comic" OR keyword.keyword = "tv-special" OR keyword.keyword = "fight" OR keyword.keyword = "violence") AND name.name LIKE "%Downey%Robert%" AND title.production_year > 2014 AND keyword.id = movie_keyword.keyword_id AND title.id = movie_keyword.movie_id AND title.id = cast_info.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND name.id = cast_info.person_id;

-- Query 99
-- 6c
SELECT keyword.keyword, name.name, title.title FROM cast_info, keyword, movie_keyword, name, title WHERE keyword.keyword = "marvel-cinematic-universe" AND name.name LIKE "%Downey%Robert%" AND title.production_year > 2014 AND keyword.id = movie_keyword.keyword_id AND title.id = movie_keyword.movie_id AND title.id = cast_info.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND name.id = cast_info.person_id;

-- Query 100
-- 6d
SELECT keyword.keyword, name.name, title.title FROM cast_info, keyword, movie_keyword, name, title WHERE (keyword.keyword = "superhero" OR keyword.keyword = "sequel" OR keyword.keyword = "second-part" OR keyword.keyword = "marvel-comics" OR keyword.keyword = "based-on-comic" OR keyword.keyword = "tv-special" OR keyword.keyword = "fight" OR keyword.keyword = "violence") AND name.name LIKE "%Downey%Robert%" AND title.production_year > 2000 AND keyword.id = movie_keyword.keyword_id AND title.id = movie_keyword.movie_id AND title.id = cast_info.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND name.id = cast_info.person_id;

-- Query 101
-- 6e
SELECT keyword.keyword, name.name, title.title FROM cast_info, keyword, movie_keyword, name, title WHERE keyword.keyword = "marvel-cinematic-universe" AND name.name LIKE "%Downey%Robert%" AND title.production_year > 2000 AND keyword.id = movie_keyword.keyword_id AND title.id = movie_keyword.movie_id AND title.id = cast_info.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND name.id = cast_info.person_id;

-- Query 102
-- 6f
SELECT keyword.keyword, name.name, title.title FROM cast_info, keyword, movie_keyword, name, title WHERE (keyword.keyword = "superhero" OR keyword.keyword = "sequel" OR keyword.keyword = "second-part" OR keyword.keyword = "marvel-comics" OR keyword.keyword = "based-on-comic" OR keyword.keyword = "tv-special" OR keyword.keyword = "fight" OR keyword.keyword = "violence") AND title.production_year > 2000 AND keyword.id = movie_keyword.keyword_id AND title.id = movie_keyword.movie_id AND title.id = cast_info.movie_id AND cast_info.movie_id = movie_keyword.movie_id AND name.id = cast_info.person_id;

-- Query 103
-- 7a
SELECT name.name, title.title FROM aka_name, cast_info, info_type, link_type, movie_link, name, person_info, title WHERE aka_name.name LIKE "%a%" AND info_type.info ="mini biography" AND link_type.link ="features" AND name.name_pcode_cf >= "A" AND name.name_pcode_cf <= "F" AND (name.gender="m" OR (name.gender = "f" AND name.name LIKE "B%")) AND person_info.note ="Volker Boehm" AND title.production_year >= 1980 AND title.production_year <= 1995 AND name.id = aka_name.person_id AND name.id = person_info.person_id AND cast_info.person_id = name.id AND title.id = cast_info.movie_id AND movie_link.linked_movie_id = title.id AND link_type.id = movie_link.link_type_id AND info_type.id = person_info.info_type_id AND person_info.person_id = aka_name.person_id AND person_info.person_id = cast_info.person_id AND aka_name.person_id = cast_info.person_id AND cast_info.movie_id = movie_link.linked_movie_id;

-- Query 104
-- 7b
SELECT name.name, title.title FROM aka_name, cast_info, info_type, link_type, movie_link, name, person_info, title WHERE aka_name.name LIKE "%a%" AND info_type.info ="mini biography" AND link_type.link ="features" AND name.name_pcode_cf LIKE "D%" AND name.gender="m" AND person_info.note ="Volker Boehm" AND title.production_year >= 1980 AND title.production_year <= 1984 AND name.id = aka_name.person_id AND name.id = person_info.person_id AND cast_info.person_id = name.id AND title.id = cast_info.movie_id AND movie_link.linked_movie_id = title.id AND link_type.id = movie_link.link_type_id AND info_type.id = person_info.info_type_id AND person_info.person_id = aka_name.person_id AND person_info.person_id = cast_info.person_id AND aka_name.person_id = cast_info.person_id AND cast_info.movie_id = movie_link.linked_movie_id;

-- Query 105
-- 7c
SELECT name.name, person_info.info FROM aka_name, cast_info, info_type, link_type, movie_link, name, person_info, title WHERE (aka_name.name LIKE "%a%" OR aka_name.name LIKE "A%") AND info_type.info ="mini biography" AND (link_type.link = "references" OR link_type.link = "referenced in" OR link_type.link = "features" OR link_type.link = "featured in") AND name.name_pcode_cf >= "A" AND name.name_pcode_cf <= "F" AND (name.gender="m" OR (name.gender = "f" AND name.name LIKE "A%")) AND title.production_year >= 1980 AND title.production_year <= 2010 AND name.id = aka_name.person_id AND name.id = person_info.person_id AND cast_info.person_id = name.id AND title.id = cast_info.movie_id AND movie_link.linked_movie_id = title.id AND link_type.id = movie_link.link_type_id AND info_type.id = person_info.info_type_id AND person_info.person_id = aka_name.person_id AND person_info.person_id = cast_info.person_id AND aka_name.person_id = cast_info.person_id AND cast_info.movie_id = movie_link.linked_movie_id;

-- Query 106
-- 8a
SELECT aka_name.name, title.title FROM aka_name, cast_info, company_name, movie_companies, name, role_type, title WHERE cast_info.note ="(voice: English version)" AND company_name.country_code ="[jp]" AND movie_companies.note LIKE "%(Japan)%" AND name.name LIKE "%Yo%" AND role_type.role ="actress" AND aka_name.person_id = name.id AND name.id = cast_info.person_id AND cast_info.movie_id = title.id AND title.id = movie_companies.movie_id AND movie_companies.company_id = company_name.id AND cast_info.role_id = role_type.id AND aka_name.person_id = cast_info.person_id AND cast_info.movie_id = movie_companies.movie_id;

-- Query 107
-- 8b
SELECT aka_name.name, title.title FROM aka_name, cast_info, company_name, movie_companies, name, role_type, title WHERE cast_info.note ="(voice: English version)" AND company_name.country_code ="[jp]" AND movie_companies.note LIKE "%(Japan)%" AND (movie_companies.note LIKE "%(2006)%" OR movie_companies.note LIKE "%(2007)%") AND name.name LIKE "%Yo%" AND role_type.role ="actress" AND title.production_year >= 2006 AND title.production_year <= 2007 AND (title.title LIKE "One Piece%" OR title.title LIKE "Dragon Ball Z%") AND aka_name.person_id = name.id AND name.id = cast_info.person_id AND cast_info.movie_id = title.id AND title.id = movie_companies.movie_id AND movie_companies.company_id = company_name.id AND cast_info.role_id = role_type.id AND aka_name.person_id = cast_info.person_id AND cast_info.movie_id = movie_companies.movie_id;

-- Query 108
-- 8c
SELECT aka_name.name, title.title FROM aka_name, cast_info, company_name, movie_companies, name, role_type, title WHERE company_name.country_code ="[us]" AND role_type.role ="writer" AND aka_name.person_id = name.id AND name.id = cast_info.person_id AND cast_info.movie_id = title.id AND title.id = movie_companies.movie_id AND movie_companies.company_id = company_name.id AND cast_info.role_id = role_type.id AND aka_name.person_id = cast_info.person_id AND cast_info.movie_id = movie_companies.movie_id;

-- Query 109
-- 8d
SELECT aka_name.name, title.title FROM aka_name, cast_info, company_name, movie_companies, name, role_type, title WHERE company_name.country_code ="[us]" AND role_type.role ="costume designer" AND aka_name.person_id = name.id AND name.id = cast_info.person_id AND cast_info.movie_id = title.id AND title.id = movie_companies.movie_id AND movie_companies.company_id = company_name.id AND cast_info.role_id = role_type.id AND aka_name.person_id = cast_info.person_id AND cast_info.movie_id = movie_companies.movie_id;

-- Query 110
-- 9a
SELECT aka_name.name, char_name.name, title.title FROM aka_name, char_name, cast_info, company_name, movie_companies, name, role_type, title WHERE (cast_info.note = "(voice)" OR cast_info.note = "(voice: Japanese version)" OR cast_info.note = "(voice) (uncredited)" OR cast_info.note = "(voice: English version)") AND company_name.country_code ="[us]" AND (movie_companies.note LIKE "%(USA)%" OR movie_companies.note LIKE "%(worldwide)%") AND name.gender ="f" AND name.name LIKE "%Ang%" AND role_type.role ="actress" AND title.production_year >= 2005 AND title.production_year <= 2015 AND cast_info.movie_id = title.id AND title.id = movie_companies.movie_id AND cast_info.movie_id = movie_companies.movie_id AND movie_companies.company_id = company_name.id AND cast_info.role_id = role_type.id AND name.id = cast_info.person_id AND char_name.id = cast_info.person_role_id AND aka_name.person_id = name.id AND aka_name.person_id = cast_info.person_id;

-- Query 111
-- 9b
SELECT aka_name.name, char_name.name, name.name, title.title FROM aka_name, char_name, cast_info, company_name, movie_companies, name, role_type, title WHERE cast_info.note = "(voice)" AND company_name.country_code ="[us]" AND movie_companies.note LIKE "%(200%)%" AND (movie_companies.note LIKE "%(USA)%" OR movie_companies.note LIKE "%(worldwide)%") AND name.gender ="f" AND name.name LIKE "%Angel%" AND role_type.role ="actress" AND title.production_year >= 2007 AND title.production_year <= 2010 AND cast_info.movie_id = title.id AND title.id = movie_companies.movie_id AND cast_info.movie_id = movie_companies.movie_id AND movie_companies.company_id = company_name.id AND cast_info.role_id = role_type.id AND name.id = cast_info.person_id AND char_name.id = cast_info.person_role_id AND aka_name.person_id = name.id AND aka_name.person_id = cast_info.person_id;

-- Query 112
-- 9c
SELECT aka_name.name, char_name.name, name.name, title.title FROM aka_name, char_name, cast_info, company_name, movie_companies, name, role_type, title WHERE (cast_info.note = "(voice)" OR cast_info.note = "(voice: Japanese version)" OR cast_info.note = "(voice) (uncredited)" OR cast_info.note = "(voice: English version)") AND company_name.country_code ="[us]" AND name.gender ="f" AND name.name LIKE "%An%" AND role_type.role ="actress" AND cast_info.movie_id = title.id AND title.id = movie_companies.movie_id AND cast_info.movie_id = movie_companies.movie_id AND movie_companies.company_id = company_name.id AND cast_info.role_id = role_type.id AND name.id = cast_info.person_id AND char_name.id = cast_info.person_role_id AND aka_name.person_id = name.id AND aka_name.person_id = cast_info.person_id;

-- Query 113
-- 9d
SELECT aka_name.name, char_name.name, name.name, title.title FROM aka_name, char_name, cast_info, company_name, movie_companies, name, role_type, title WHERE (cast_info.note = "(voice)" OR cast_info.note = "(voice: Japanese version)" OR cast_info.note = "(voice) (uncredited)" OR cast_info.note = "(voice: English version)") AND company_name.country_code ="[us]" AND name.gender ="f" AND role_type.role ="actress" AND cast_info.movie_id = title.id AND title.id = movie_companies.movie_id AND cast_info.movie_id = movie_companies.movie_id AND movie_companies.company_id = company_name.id AND cast_info.role_id = role_type.id AND name.id = cast_info.person_id AND char_name.id = cast_info.person_role_id AND aka_name.person_id = name.id AND aka_name.person_id = cast_info.person_id;

