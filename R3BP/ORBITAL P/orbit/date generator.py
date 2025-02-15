# date format required    2019-02-19T18:00:00

from datetime import datetime, timedelta
import random

# Define the start date and time interval
start_date = datetime(2022, 8, 10, 0, 0, 0)  # Start date and time
time_interval = timedelta(days=60)  # Time interval of 30 days (adjust as needed)

# Number of dates to generate
num_dates = 100 # You can change this number as needed

# Generate random dates
random_dates = [start_date + i * time_interval for i in range(num_dates)]

# Format the dates as strings with apostrophes and commas
formatted_dates = ["'" + date.strftime('%Y-%m-%d %H:%M:%S') + "'," for date in random_dates]

# Print the list of random dates
for date_str in formatted_dates:
    print(date_str)
