/*
 * IBMisc: Misc. Routines for Ibmisc (and other code)
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <ibmisc/array.hpp>

namespace ibmisc {

// =============================================================
class Datetime;
class Time;

struct Date : public std::array<int,3>
{
    typedef std::array<int,3> super;

    int &year() { return (*this)[0]; }
    int &month() { return (*this)[1]; }
    int &day() { return (*this)[2]; }

    int year() const { return (*this)[0]; }
    int month() const { return (*this)[1]; }
    int day() const { return (*this)[2]; }

    Date(int _year=0, int _month=1, int _day=1) :
        super(make_array(_year, _month, _day)) {}

    Date(Datetime const &ymdt);

    Datetime operator+(Time const &tm);

};

struct Time : public std::array<int,4>
{
    typedef std::array<int,4> super;

    int &hour() { return (*this)[0]; }
    int &minute() { return (*this)[1]; }
    int &second() { return (*this)[2]; }
    int &usec() { return (*this)[3]; }

    int hour() const { return (*this)[0]; }
    int minute() const { return (*this)[1]; }
    int second() const { return (*this)[2]; }
    int usec() const { return (*this)[3]; }


    Time(int _hour=0, int _minute=0, int _second=0, int _usec=0)
        : super(make_array(_hour, _minute, _second, _usec)) {}

    Time(Datetime const &ymdt);
};

struct Datetime : public std::array<int,7>
{
    typedef std::array<int,7> super;

    int &year() { return (*this)[0]; }
    int &month() { return (*this)[1]; }
    int &day() { return (*this)[2]; }
    int &hour() { return (*this)[3]; }
    int &minute() { return (*this)[4]; }
    int &second() { return (*this)[5]; }
    int &usec() { return (*this)[6]; }

    int year() const { return (*this)[0]; }
    int month() const { return (*this)[1]; }
    int day() const { return (*this)[2]; }
    int hour() const { return (*this)[3]; }
    int minute() const { return (*this)[4]; }
    int second() const { return (*this)[5]; }
    int usec() const { return (*this)[6]; }

//    Datetime() : Datetime(0,1,1,0,0,0,0) {}

    Datetime(
        int _year=0, int _month=1, int _day=1,
        int _hour=0, int _minute=0, int _second=0, int _usec=0) :
        super(make_array(_year, _month, _day, _hour, _minute, _second, _usec)) {}

    Datetime(Date const &ymd, Time const &tm) {
        for (int i=0; i<3; ++i) (*this)[i] = ymd[i];
        for (int i=0; i<4; ++i) (*this)[i+3] = tm[i];
    }

    Datetime(Date const &ymd) {
        for (int i=0; i<3; ++i) (*this)[i] = ymd[i];
        for (int i=0; i<4; ++i) (*this)[i+3] = 0;
    }

    Datetime(Time const &tm) {
        for (int i=0; i<3; ++i) (*this)[i] = 0;
        for (int i=0; i<4; ++i) (*this)[i+3] = tm[i];
    }

};


inline Date::Date(Datetime const &ymdt)
    { for (int i=0; i<3; ++i) (*this)[i] = ymdt[i]; }

inline Time::Time(Datetime const &ymdt)
    { for (int i=0; i<4; ++i) (*this)[i] = ymdt[i+3]; }

inline Datetime Date::operator+(Time const &tm)
    { return Datetime(*this, tm); }

// =============================================================
// ---------------------------------------------------
/** Julian day or whatnot */
struct JDate : public std::array<long,1>
{
    typedef std::array<long,1> super;
    long &day() { return (*this)[0]; }
    long day() const { return (*this)[0]; }

    JDate(long _val)
        : super(make_array(_val)) {}

    JDate operator+(long nday) const { return JDate(day() + nday); }
    JDate operator-(long nday) const { return JDate(day() - nday); }
};

/** usec within a day */
struct JTime : public std::array<long,1>
{
    typedef std::array<long,1> super;
    long &usec() { return (*this)[0]; }
    long usec() const { return (*this)[0]; }

    JTime(long _val) : super(make_array(_val)) {}
};

struct JDatetime : public std::tuple<JDate, JTime> {
    typedef std::tuple<JDate, JTime> super;
    JDate &jdate() { return std::get<0>(*this); }
    JTime &jtime() { return std::get<1>(*this); }

    JDate jdate() const { return std::get<0>(*this); }
    JTime jtime() const { return std::get<1>(*this); }

    JDatetime(JDate _jdate, JTime _jtime) :
        super(std::make_tuple(_jdate, _jtime)) {}
};

// -------------------------------------------------
// JTime Conversions
struct USEC {
    static const long DAY = 86400L*1000000L;
    static const long HOUR = 3600L*1000000L;
    static const long MINUTE = 60L*1000000L;
    static const long SECOND = 1000000L;
};

inline JTime to_jtime(Time const &t) {
    return JTime(
        t.hour()*USEC::HOUR +
        t.minute()*USEC::MINUTE +
        t.second()*USEC::SECOND +
        t.usec());
}
inline Time to_time(JTime jt) {
    return Time(
        jt.usec() / USEC::HOUR,
        (jt.usec() % USEC::HOUR) / USEC::MINUTE,
        (jt.usec() % USEC::MINUTE) / USEC::SECOND,
        jt.usec() % USEC::SECOND);
}

// ---------------------------------------------------
// Calendar Conversions

struct Calendar {
    virtual ~Calendar();

    /** Used in the NetCDF `calendar` attribute (eg: "365_day"). */
    virtual std::string to_cf() const = 0;

    virtual JDate to_jdate(Date const &da) const = 0;
    virtual Date to_date(JDate jda) const = 0;

    JDatetime to_jdatetime(Datetime const &dt) const
        { return JDatetime(to_jdate(dt), to_jtime(dt)); }

    Datetime to_datetime(JDatetime const &jdt) const
        { return Datetime(to_date(jdt.jdate()), to_time(jdt.jtime())); }

    Date add(Date const da, long nday) const
        { return to_date(to_jdate(da) + nday); }
    Date sub(Date const da, long nday) const
        { return to_date(to_jdate(da) - nday); }

};

struct Cal365 : public Calendar {
    static const std::array<std::string,12> month_name3;
    static const std::array<int,12> month_len;
    static const std::array<int,13> month_start;

    virtual ~Cal365();

    std::string to_cf() const;
    JDate to_jdate(Date const &cda) const;
    Date to_date(JDate da) const;
};
// ---------------------------------------------------
extern Cal365 cal365;

/** Encodes units like "seconds since 1950-01-01 15:30:00" */
struct TimeUnit {
    enum TimeElement {YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, USEC};

    Calendar const *cal;    // Calendar we use
    Datetime base;
    TimeElement uniti;    // 0=years, 6=usec
    JDatetime basej;    // Julian day and seconds-in-day

public:
    TimeUnit();

    TimeUnit(Calendar const *_cal, Datetime const &_base, TimeElement _uniti);

    /** Produces a CF-compliant string describing this unit for time. */
    std::string to_cf() const;

    JDatetime to_jdatetime(double tm);

    Datetime to_datetime(double tm);

};



}    // namespace

extern std::ostream &operator<<(std::ostream &out, ibmisc::Date const &da);
extern std::ostream &operator<<(std::ostream &out, ibmisc::Time const &da);
extern std::ostream &operator<<(std::ostream &out, ibmisc::Datetime const &da);

