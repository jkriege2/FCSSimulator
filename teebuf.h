/*
    Copyright (c) 2008-2015 Jan W. Krieger (<jan@jkrieger.de>), German Cancer Research Center + IWR, University Heidelberg

    This software is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef TEEBUF_H
#define TEEBUF_H



#include <cstdio>
#include <streambuf>
#include <iostream>

// source: http://wordaligned.org/articles/cpp-streambufs
class teebuf: public std::streambuf
{
public:
    // Construct a streambuf which tees output to both input
    // streambufs.
    teebuf(std::streambuf * sb1, std::streambuf * sb2)
        : sb1(sb1)
        , sb2(sb2)
    {
    }
    std::streambuf * get_sb1() const { return sb1; };
    std::streambuf * get_sb2() const { return sb2; };
private:
    // This tee buffer has no buffer. So every character "overflows"
    // and can be put directly into the teed buffers.
    virtual int overflow(int c)
    {
        if (c == EOF)
        {
            return !EOF;
        }
        else
        {
            int const r1 = sb1->sputc(c);
            int const r2 = sb2->sputc(c);
            return r1 == EOF || r2 == EOF ? EOF : c;
        }
    }

    // Sync both teed buffers.
    virtual int sync()
    {
        int const r1 = sb1->pubsync();
        int const r2 = sb2->pubsync();
        return r1 == 0 && r2 == 0 ? 0 : -1;
    }

private:
    std::streambuf * sb1;
    std::streambuf * sb2;
};


class ScopedDelayedStreamDoubler2 {
public:
    ScopedDelayedStreamDoubler2()
    { redirected=false; mOriginal=NULL; mOldBuffer=NULL; tee=NULL; }

    void redirect(std::ostream & inOriginal, std::ostream & inRedirect) {
        if (redirected) {
            mOriginal->rdbuf(tee->get_sb1());
            delete tee;
            tee=NULL;
            redirected=false;
        }
        redirected=true;
        tee=new teebuf(inOriginal.rdbuf(), inRedirect.rdbuf());
        mOldBuffer=inOriginal.rdbuf(tee);
        mOriginal=&inOriginal;
    }

    ~ScopedDelayedStreamDoubler2()
    {
        if (redirected)  {
            mOriginal->rdbuf(tee->get_sb1());
            delete tee;
            tee=NULL;
            redirected=false;
        }
    }

private:
    ScopedDelayedStreamDoubler2(const ScopedDelayedStreamDoubler2&);
    ScopedDelayedStreamDoubler2& operator=(const ScopedDelayedStreamDoubler2&);

    bool redirected;
    std::ostream* mOriginal;
    std::streambuf * mOldBuffer;
    teebuf* tee;
};

#endif // TEEBUF_H
