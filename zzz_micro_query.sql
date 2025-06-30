SELECT table_1.id, table_1.col_1, table_2.id, table_2.col_2 FROM table_1, table_2 WHERE table_1.id = table_2.id;


SELECT table_3.id, table_3.col_3, table_1.id, table_1.col_1 FROM table_3, table_1 WHERE table_3.id = table_1.id;
SELECT table_2.id, table_2.col_2, table_1.id, table_1.col_1 FROM table_2,  table_1 WHERE table_2.id = table_1.id ;
SELECT table_1.id, table_1.col_1, table_2.id, table_2.col_2, table_3.id, table_3.col_3 FROM table_1, table_2, table_3 WHERE table_1.id = table_2.id AND table_1.id = table_3.id;
