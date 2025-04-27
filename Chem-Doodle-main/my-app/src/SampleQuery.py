from cassandra.cluster import Cluster

cluster = Cluster(['96.125.118.113'], port=9042)
session = cluster.connect('your_keyspace')

rows = session.execute('SELECT * FROM your_table LIMIT 10')
for row in rows:
    print(row)
