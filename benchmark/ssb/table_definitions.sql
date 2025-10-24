CREATE DATABASE ssb;
USE ssb;

-- dimension: customer
CREATE TABLE customer (
    c_custkey     INT(4)      NOT NULL,
    c_name        VARCHAR(25) NOT NULL,
    c_address     VARCHAR(25) NOT NULL,
    c_city        VARCHAR(10) NOT NULL,
    c_nation      VARCHAR(15) NOT NULL,
    c_region      VARCHAR(12) NOT NULL,
    c_phone       VARCHAR(15) NOT NULL,
    c_mktsegment  VARCHAR(10) NOT NULL
);

-- dimension: supplier
CREATE TABLE supplier (
    s_suppkey     INT(4)      NOT NULL,
    s_name        VARCHAR(25) NOT NULL,
    s_address     VARCHAR(25) NOT NULL,
    s_city        VARCHAR(10) NOT NULL,
    s_nation      VARCHAR(15) NOT NULL,
    s_region      VARCHAR(12) NOT NULL,
    s_phone       VARCHAR(15) NOT NULL
);

-- dimension: part
CREATE TABLE part (
    p_partkey     INT(4)      NOT NULL,
    p_name        VARCHAR(22) NOT NULL,
    p_mfgr        CHAR(6)     NOT NULL,
    p_category    CHAR(7)     NOT NULL,
    p_brand1      CHAR(9)     NOT NULL,
    p_color       VARCHAR(11) NOT NULL,
    p_type        VARCHAR(25) NOT NULL,
    p_size        INT(4)      NOT NULL,
    p_container   CHAR(10)    NOT NULL
);

-- dimension: dwdate (date dimension)
CREATE TABLE dwdate (
    d_datekey          INT(4)   NOT NULL,
    d_date             DATE     NOT NULL,
    d_dayofweek        VARCHAR(9),
    d_month            VARCHAR(9),
    d_year             INT(4),
    d_yearmonthnum     INT(4),
    d_yearmonth        VARCHAR(7),
    d_daynuminweek     INT(4),
    d_daynuminmonth    INT(4),
    d_daynuminyear     INT(4),
    d_monthnuminyear   INT(4),
    d_weeknuminyear    INT(4),
    d_sellingseason    VARCHAR(13),
    d_lastdayinmonthfl CHAR(1),
    d_holidayfl        CHAR(1),
    d_weekdayfl        CHAR(1)
);

-- fact: lineorder
CREATE TABLE lineorder (
    lo_orderkey       INT(4)   NOT NULL,
    lo_linenumber     INT(4)   NOT NULL,
    lo_custkey        INT(4)   NOT NULL,
    lo_partkey        INT(4)   NOT NULL,
    lo_suppkey        INT(4)   NOT NULL,
    lo_orderdate      DATE     NOT NULL,
    lo_orderpriority  VARCHAR(15) NOT NULL,
    lo_shippriority   CHAR(1)  NOT NULL,
    lo_quantity       INT(4)   NOT NULL,
    lo_extendedprice  INT(4)   NOT NULL,
    lo_ordtotalprice  INT(4)   NOT NULL,
    lo_discount       INT(4)   NOT NULL,
    lo_revenue        INT(4)   NOT NULL,
    lo_supplycost     INT(4)   NOT NULL,
    lo_tax            INT(4)   NOT NULL,
    lo_commitdate     DATE     NOT NULL,
    lo_shipmode       VARCHAR(10) NOT NULL
);
