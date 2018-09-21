#include<string>
#include<time.h>
#include<math.h>
#include<iostream>
#include "ConstDefine.h"

// very basic functions

using namespace std;
bool IsLeap(int year)
{
	return (year % 4 == 0 || year % 400 == 0) && (year % 100 != 0);// wrong ͬһ�����ԶԽ����Ӱ��
	// return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}

//�����StringToDate��������ȡ�������е������ղ��ж������Ƿ�Ϸ�
//���ַ������������ �涨���ڵĸ�ʽ��yyyy-mm-dd
bool StringToDate(string date, int& year, int& month, int& day)
{
	year = atoi((date.substr(0, 4)).c_str());
	month = atoi((date.substr(5, 2)).c_str());
	day = atoi((date.substr(8, 2)).c_str());
	int DAY[12] = { 31,28,31,30,31,30,31,31,30,31,30,31 };
	if (IsLeap(year)) {
		DAY[1] = 29;
	}
	return year >= 0 && month <= 12 && month>0 && day <= DAY[month - 1] && day>0;
}

//DayInYear�ܸ��ݸ��������ڣ�������ڸ���ĵڼ��죬��������
int DayInYear(int year, int month, int day)
{
	//int _day = 0;
	int DAY[12] = { 31,28,31,30,31,30,31,31,30,31,30,31 };
	if (IsLeap(year))
		DAY[1] = 29;
	if (month == 1)
		return day;
	else
	{
		for (int i = 0; i < month - 1; ++i)
		{
			day += DAY[i];
		}
	}
	return day;
}
int DaysBetween2Date(string date1, string date2)
{
	//ȡ�������е�������
	int year1, month1, day1;
	int year2, month2, day2;
	if (!StringToDate(date1, year1, month1, day1) || !StringToDate(date2, year2, month2, day2))
	{
		cout << "��������ڸ�ʽ����ȷ";
		return -1;
	}

		//���¾���ͬ

	if (year1 == year2 && month1 == month2)
	{
		return day1 > day2 ? day1 - day2 : day2 - day1;

		//����ͬ
	}
	else if (year1 == year2)
	{
		int d1, d2;
		d1 = DayInYear(year1, month1, day1);
		d2 = DayInYear(year2, month2, day2);
		return d1 > d2 ? d1 - d2 : d2 - d1;

		//���¶�����ͬ
	}
	else {
		//ȷ��year1��ݱ�year2��
		if (year1 > year2)
		{
			//swap��������ֵ�Ľ���
			swap(year1, year2);
			swap(month1, month2);
			swap(day1, day2);
		}
		int d1, d2, d3;
		if (IsLeap(year1))
			d1 = 366 - DayInYear(year1, month1, day1); //ȡ����������ڸ��껹���¶�����
		else
			d1 = 365 - DayInYear(year1, month1, day1);
		d2 = DayInYear(year2, month2, day2); //ȡ���ڵ����еĵڼ���
											 // cout << "d1:" << d1 << ", d2:" << d2;

		d3 = 0;
		for (int year = year1 + 1; year < year2; year++)
		{
			if (IsLeap(year))
				d3 += 366;
			else
				d3 += 365;
		}
		return d1 + d2 + d3;
	}
}
