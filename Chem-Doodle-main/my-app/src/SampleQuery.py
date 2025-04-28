from cassandra.cluster import Cluster

cluster = Cluster(['96.125.118.113'], port=9042)
session = cluster.connect('drug_discovery')

rows = session.execute('SELECT * FROM compounds LIMIT 10')
for row in rows:
    print(row)
