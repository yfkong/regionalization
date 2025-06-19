Areal data file format：
-----------------------------------------
186	3	4	1
1	1	1	1	2	1	3
id	x	y	Area	attr1	attr2	attr2
1	4.2	1.7	0.97	0.78	0.58	0.33
......
186	0.2	3.4	0.76	0.36	0.45	0.14
-----------------------------------------
1st line: number of spatial units (186), number of attributes (3), begin index of attribtes (4), and data standardlized (1 for yes);
2nd line: weights of attrbues, 2 on attr1, 1 on attr2, and 3 on attr3;
3rd line: field names of data: id, x, y, area, and attribute names.
4th and following lines: data.

Connectivity data file format：
-----------------------------------
Idx	id1	id2
1	1	5
2	1	6
3	2	6
4	2	8
...
997	186	180
-----------------------------------
In each line, unit id1 and  unit id2 are neighbors.

Run the Python program:
\pypy2\pypy rg2025.py g120_5a.txt g120_conn.txt 0 0 5 1 10
