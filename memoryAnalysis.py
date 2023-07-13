import tracemalloc
tracemalloc.start(10)

from subprocess import call
call(["python", "main.py"])

snapshot2 = tracemalloc.take_snapshot()
#top_stats = snapshot.statistics('lineno')

top_stats = snapshot2.compare_to(snapshot1, 'lineno')

print("[ Top 10 ]")
for stat in top_stats[:10]:
    print(stat)
