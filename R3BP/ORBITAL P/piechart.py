import matplotlib.pyplot as plt

# Names of the major contributors
labels = ['Ceres', 'Vesta', 'Pallas', 'Hygiea', 'Other']

# Mass distribution percentages
# Ceres: 31%, Vesta+Pallas+Hygiea: 20%, Others: 49%
mass_distribution = [31,9,7,3, 49]

# Colors for different sections
colors = ['#ff9999','#66b3ff','#FFA07A','#add8e6','#98fb98']

# Exploding the 1st slice (Ceres)
explode = (0,0,0,0,0)

# Plotting the pie chart
plt.figure(figsize=(8, 8))
plt.pie(mass_distribution, labels=labels, colors=colors, explode=explode, autopct='%1.1f%%', startangle=140)

#plt.title('Mass Distribution of the Main Asteroid Belt')
plt.axis('equal')  # Equal aspect ratio ensures the pie chart is circular.

# Display the chart
plt.show()
