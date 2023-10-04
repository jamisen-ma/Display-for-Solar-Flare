import bisect
import math
from datetime import datetime


def convert_timedelta(duration):
    seconds = duration.seconds
    minutes = math.ceil(seconds / 60) + 1
    return minutes

def to_datetime_object(date_string, date_format):
    s = datetime.strptime(str(date_string), date_format)
    return s


def convert_timedelta(duration):
    seconds = duration.seconds
    minutes = ((seconds % 3600) // 60) + 1
    return minutes


def roundup(x):
    return int(math.ceil(x / 100.0)) * 100


def rounddown(x):
    return (int(math.ceil(x / 100.0)) * 100) - 100