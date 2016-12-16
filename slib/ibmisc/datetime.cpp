#include <cmath>
#include <array>
#include <algorithm>
#include <boost/format.hpp>
#include <ibmisc/datetime.hpp>

namespace ibmisc {

// ===================================================
// Cal365
Calendar::~Calendar() {}

static std::array<int,13> init_month_start()
{
    std::array<int,13> ret;
    int start = 0;
    for (int i=0; i<13; ) {
        start += Cal365::month_len[i];
        ++i;
        ret[i] = start;
    }
    return ret;
}
const std::array<std::string,12> Cal365::month_name3
    {"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};
const std::array<int,12> Cal365::month_len {31,28,31,30,31,30,31,31,30,31,30,31};
const std::array<int,13> Cal365::month_start = init_month_start();

Cal365 cal365;    // Singleton

// ------------------------------------------------------
/** Produces a Juliay Day-like number. */
JDate Cal365::to_jdate(Date const &da) const
{
    return JDate(da.year()*365L + month_start[da.month()-1] + (da.day()-1));
}

Date Cal365::to_date(JDate jda) const
{
    int yy = jda.day() / 365;

    // Returns an iterator pointing to the first element in the range
    // [first, last) that is greater than value.
    int day = jda.day() % 365;
    auto mm(std::upper_bound(month_start.begin(), month_start.end(), day) - 1);
    int dd = day - *mm;

    return Date(yy,mm-month_start.begin()+1,dd+1);
}

// ------------------------------------------------------

// =====================================================
// TimeUnit
TimeUnit::TimeUnit(Calendar const *_cal, Datetime const &_base, TimeUnit::TimeElement _uniti)
    : cal(_cal), base(_base), uniti(_uniti), basej(cal->to_jdatetime(base))
{}

TimeUnit::TimeUnit()
    : cal(&cal365), base(1970,1,1), uniti(TimeElement::SECOND),
    basej(cal->to_jdatetime(base))
{}


JDatetime TimeUnit::to_jdatetime(double tm)
{
    switch(uniti) {
    case TimeElement::SECOND : {
        // Split floating point into integer and fractional parts
        double tm_int0, tm_frac;
        tm_frac = modf(tm, &tm_int0);
        long tm_int = tm_int0;

        // Determine days and usec-with-day for the integer part
        long tm_int_day = tm_int / 86400L;
        long tm_int_usec = (tm_int % 86400L) * 1000000L;

        // Determine extra usec for the fractional part
        long tm_frac_usec = (tm_frac * 1000000. + .5);

        // Determine total usec from all portions
        long usec = basej.jtime().usec() + tm_int_usec + tm_frac_usec;

        return JDatetime(
            basej.jdate() + tm_int_day + (usec / USEC::DAY),
            JTime(usec % USEC::DAY));
    }
    }
}

std::string TimeUnit::to_cf()
{
    switch(uniti) {
    case TimeElement::SECOND : {
        return (boost::format("seconds since %04d-%02d-%02d %02d:%02d:%02d")
            % base.year() % base.month() % base.day()
            % base.hour() % base.minute() % base.second()).str();
    }}
}


Datetime TimeUnit::to_datetime(double tm)
{
    return cal->to_datetime(to_jdatetime(tm));
}


}    // namespace

using namespace ibmisc;

std::ostream &operator<<(std::ostream &os, ibmisc::Date const &da)
{
    os << (boost::format
        ("%04d-%02d-%02d") % da.year() % da.month() % da.day()).str();
    return os;
}

std::ostream &operator<<(std::ostream &os, ibmisc::Time const &da)
{
    os << (boost::format("%02d:%02d") % da.hour() % da.minute());
    if (da.second() != 0 || da.usec() != 0) os << (boost::format(":%02d") % da.second());
    if (da.usec() != 0)  os << (boost::format(".%06d") % da.usec());
    return os;
}
std::ostream &operator<<(std::ostream &os, ibmisc::Datetime const &dt)
{
    os << Date(dt) << "T" << Time(dt);
    return os;
}
