//
// Created by Chuang on 2019/5/15.
//

#pragma once

#define PeriodLen 50.0

#include <string>
#include <sstream>

template <class Type>
Type stringToNum(const std::string& str)
{
    std::istringstream iss(str);
    Type num;
    iss >> num;
    return num;
}
double naivetime(std::string l);
int getPeriod(double time);
int getMaxPeriod();
double getPeriodStart(double time);
double getPeriodEnd(double time);