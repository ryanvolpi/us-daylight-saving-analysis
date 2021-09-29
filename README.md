# US Daylight Saving Analysis

This analysis investigates the effect Daylight Saving Time (DST) has on the american population's sunrise and sunset times. The purpose of doing this is to support my position in arguments with my friends that DST should be extended year round.

To understand the effect of DST on the distribution of sunrise and sunset times in the United States, I estimate what time the sun rises and sets every day for each person in the country. The time differs from person to person due to changes in geographical location (eg the on 7/10/2021 the sunrises at ). We define geographical location at the level of census block. There are over 11 million census blocks and 365 days in the year, which equates to over 4 billion unique sunrise and sunset times. The primary challenge of this project, beyond identifying and joining relevant data sources, is calculating these 4+ billion sunrise/sunset times in a time-efficient manner.

The steps are the following:

1. Download shape file and population data for each census block
2. Calculate centroid (latitude / longitude) for each census block
3. Get timezone associated with each census block
4. Calculate sunrise / sunset time for each day for each census block based on centroid
5. Create plot of sunrise/sunset distributions with and without DST

![image](https://user-images.githubusercontent.com/22488709/135197798-edffced6-98f8-43bf-aa7b-414ed8f5b819.png)
