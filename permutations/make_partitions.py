def make_partitions(elements):
    yield from _make_partitions(sorted(elements, reverse=True), [], []) 

def _make_partitions(elements, active_partitions, inactive_partitions):
    if not elements:
        yield active_partitions + inactive_partitions
        return

    elem = elements.pop()

    # Make create a new partition
    active_partitions.append([elem])
    yield from _make_partitions(elements, active_partitions, inactive_partitions)
    active_partitions.pop()

    # Add element to each existing partition in turn
    size = len(active_partitions)
    for part in active_partitions[::-1]:
        part.append(elem)
        yield from _make_partitions(elements, active_partitions, inactive_partitions)
        part.pop()

        # Remove partition that would create a cross if new elements were added
        inactive_partitions.append(active_partitions.pop())
    
    # Add back removed partitions
    for _ in range(size):
        active_partitions.append(inactive_partitions.pop())

    elements.append(elem)

for partitions in make_partitions([1, 2, 3, 4, 5]):
    print(partitions[0])