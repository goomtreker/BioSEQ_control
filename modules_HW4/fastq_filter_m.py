def mean(array: list[int | float]) -> float:
    return sum(array)/len(array)


def phread_score(field_4: str) -> float:  # Функция для расчета phread score
    symbol_collector = []
    for symbol in field_4:
        symbol_collector.append(ord(symbol) - 33)
    return mean(symbol_collector)


def check_limits(gc_bounds: tuple | int, length_bounds: tuple | int) -> tuple:
    """Определение лимитов,
    в том случае если подали
    одним числом"""
    if type(gc_bounds) is int:
        gc_bounds = (0, gc_bounds)
    elif type(length_bounds) is int:
        length_bounds = (0, length_bounds)
    return gc_bounds, length_bounds
