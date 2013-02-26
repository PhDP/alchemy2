/***************************************************************************
*   Copyright (C) 2004 by Vibhav Gogate                                   *
*   vgogate@ics.uci.edu                                                   *
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
*   This program is distributed in the hope that it will be useful,       *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
*   GNU General Public License for more details.                          *
*                                                                         *
*   You should have received a copy of the GNU General Public License     *
*   along with this program; if not, write to the                         *
*   Free Software Foundation, Inc.,                                       *
*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
**************************************************************************/
#ifndef UTIL_HPP
#define UTIL_HPP
#include <iostream>
#include <vector>
#include <algorithm>


using namespace std;

template<typename _Classname>
void do_set_union(const vector<_Classname> &a, const vector<_Classname> &b, vector<_Classname> &c_)
{
	vector<_Classname> c;
	c.resize(a.size()+b.size());
	typename vector<_Classname>::iterator curr,last;
	last=c.end();
	curr=set_union(a.begin(),a.end(),b.begin(),b.end(),c.begin());
	c.erase(curr,last);
	c_=c;
}

template<typename _Classname, typename _Compare>
void do_set_union(const vector<_Classname> &a, const vector<_Classname> &b, vector<_Classname> &c_,_Compare comp)
{
	vector<_Classname> c;
	c.resize(a.size()+b.size());
	typename vector<_Classname>::iterator curr,last;
	last=c.end();
	curr=set_union(a.begin(),a.end(),b.begin(),b.end(),c.begin(),comp);
	c.erase(curr,last);
	c_=c;
}

// see if a is included in b

template<typename _Classname>
bool do_set_inclusion(const vector<_Classname> &a, const vector<_Classname> &b)
{
	return includes(b.begin(),b.end(),a.begin(),a.end());
}
template<typename _Classname>
bool do_set_inclusion(const _Classname &a_, const vector<_Classname> &b)
{
	vector<_Classname> a;
	a.push_back(a_);
	return includes(b.begin(),b.end(),a.begin(),a.end());
}

template<typename _Classname, typename _Compare>
bool do_set_inclusion(const vector<_Classname> &a, const vector<_Classname> &b, _Compare comp)
{
	return includes(b.begin(),b.end(),a.begin(),a.end(),comp);
}
template<typename _Classname, typename _Compare>
bool do_set_inclusion(const _Classname &a_, const vector<_Classname> &b,_Compare comp)
{
	vector<_Classname> a;
	a.push_back(a_);
	return includes(b.begin(),b.end(),a.begin(),a.end(),comp);
}

template<typename _Classname>
bool do_set_intersection(const vector<_Classname> &a, const vector<_Classname> &b, vector<_Classname> &c_)
{
	vector<_Classname> c;
	c.resize(a.size());
	typename vector<_Classname>::iterator curr,last;
	last=c.end();
	curr=set_intersection(a.begin(),a.end(),b.begin(),b.end(),c.begin());
	c.erase(curr,last);
	c_=c;
	return !c_.empty();
}

template<typename _Classname, typename _Compare>
bool do_set_intersection(const vector<_Classname> &a, const vector<_Classname> &b, vector<_Classname> &c_,_Compare comp)
{
	vector<_Classname> c;
	c.resize(a.size());
	typename vector<_Classname>::iterator curr,last;
	last=c.end();
	curr=set_intersection(a.begin(),a.end(),b.begin(),b.end(),c.begin(),comp);
	c.erase(curr,last);
	c_=c;
	return !c_.empty();
}

template<typename _Classname>
void do_set_difference(const vector<_Classname> &a, const vector<_Classname> &b, vector<_Classname> &c_)
{
	vector<_Classname> c;
	c.resize (a.size());
	typename vector<_Classname>::iterator curr,last;
	last=c.end();
	curr=set_difference(a.begin(),a.end(),b.begin(),b.end(),c.begin());
	c.erase(curr,last);
	c_=c;
}

template<typename _Classname, typename _Compare>
void do_set_difference(const vector<_Classname> &a, const vector<_Classname> &b, vector<_Classname> &c_,_Compare comp)
{
	vector<_Classname> c;
	c.resize(a.size());
	typename vector<_Classname>::iterator curr,last;
	last=c.end();
	curr=set_difference(a.begin(),a.end(),b.begin(),b.end(),c.begin(),comp);
	c.erase(curr,last);
	c_=c;
}

bool is_disjoint(const vector<int> &a, const vector<int> &b)
{
	if(a.empty() || b.empty()) return true;
	vector<int> c;
	do_set_intersection(a,b,c);
	return c.empty();
}
#endif
