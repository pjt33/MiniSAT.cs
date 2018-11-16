/******************************************************************************************
MiniSat -- Copyright (c) 2003-2005, Niklas Een, Niklas Sorensson
MiniSatCS -- Copyright (c) 2006-2007 Michal Moskal

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

using System;
using System.IO;
using System.Text;
using System.Diagnostics;
using System.Collections.Generic;

// NOTE! Variables are just integers. No abstraction here. They should be chosen from 0..N,
// so that they can be used as array indices.
using Var = System.Int32;

namespace MiniSatCS {

public class Solver {

	#region lbool
	// Lifted booleans:

	// the problem is C# allows use of ~ on any enum type
	// thefore ~lbool.True == lbool.False, but we also end up using 
	// two undef values
	public enum lbool : sbyte {
		True = 1,
		False = -2,
		Undef0 = 0,
		Undef1 = -1
	}

	public const lbool l_True  = lbool.True;
	public const lbool l_False = lbool.False;

	static lbool toLbool(bool v) { return v ? lbool.True : lbool.False; }

	public static bool isUndef(lbool l)
	{
		return l != lbool.True && l != lbool.False;
	}
	#endregion

	#region Literals
	//=================================================================================================
	// Variables, literals, clause IDs:

	const int var_Undef = -1;


	public struct Lit {
		public int     x;
		//TODO we cannot do that, is that a problem?
		//public Lit() : x(2*var_Undef) {}   // (lit_Undef)
		public Lit(Var var, bool sign) {
		  x = var + var + (sign ? 1 : 0);
		}
		public Lit(Var var) {
		  x = var + var;
		}
		public static Lit operator ~ (Lit p) { Lit q; q.x = p.x ^ 1; return q; }
		public static bool operator == (Lit p, Lit q) { return index(p) == index(q); }
		public static bool operator != (Lit p, Lit q) { return index(p) != index(q); }
		public static bool operator <  (Lit p, Lit q) { return index(p)  < index(q); }  // '<' guarantees that p, ~p are adjacent in the ordering.

		public override string ToString()
		{
			return (sign(this) ? "-" : "") + "x" + var(this);
		}

		public override int GetHashCode()
		{
			return x;
		}

		public override bool Equals(object other)
		{
			if (other == null) return false;

			if (other is Lit)
				return (Lit)other == this;
			return false;
		}
	}


	static  public bool sign  (Lit p) { return (p.x & 1) != 0; }
	static  public int  var   (Lit p) { return p.x >> 1; }
	static  public int  index (Lit p) { return p.x; }                // A "toInt" method that guarantees small, positive integers suitable for array indexing.
	//static  Lit  toLit (int i) { Lit p = new Lit(); p.x = i; return p; }  // Inverse of 'index()'.
	//static  Lit  unsign(Lit p) { Lit q = new Lit(); q.x = p.x & ~1; return q; }
	//static  Lit  id    (Lit p, bool sgn) { Lit q; q.x = p.x ^ (sgn ? 1 : 0); return q; }

	static public Lit lit_Undef = new Lit(var_Undef, false);  // \- Useful special constants.
	//static Lit lit_Error = new Lit(var_Undef, true );  // /
	#endregion

    #region Clauses
	//=================================================================================================
	// Clause -- a simple class for representing a clause:

	const   int ClauseId_null = int.MinValue;

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	public class Clause {
		Lit[]     data;
		bool is_learnt;
		protected internal Clause(bool learnt, vec<Lit> ps) {
			is_learnt = learnt;
			this.data = new Lit[ps.size()];
			for (int i = 0; i < ps.size(); i++) data[i] = ps[i];
		}


		public int       size        ()       { return data.Length; }
		public bool      learnt      ()       { return is_learnt; }
		
		public Lit this[int i]
		{
		  get { return data[i]; }
		  set { data[i] = value; }
		}

		public float activity;

		public override string ToString()
		{
			StringBuilder sb = new StringBuilder();
			sb.Append("[");
			foreach (Lit l in data)
				sb.Append(l).Append(", ");
			sb.Append("]");
			return sb.ToString();
		}

		public Lit[] GetData() { return data; }
	}
		
	protected virtual Clause Clause_new(bool learnt, vec<Lit> ps) 
	{
		return new Clause(learnt, ps); 
	}

	#endregion

	#region Utilities
	//=================================================================================================
	// Random numbers:

	// Returns a random float 0 <= x < 1. Seed must never be 0.
	static double drand(ref double seed) 
	{
		seed *= 1389796;
		int q = (int)(seed / 2147483647);
		seed -= (double)q * 2147483647;
		return seed / 2147483647; 
	}

	// Returns a random integer 0 <= x < size. Seed must never be 0.
	static int irand(ref double seed, int size) 
	{
		return (int)(drand(ref seed) * size);
	}

	//=================================================================================================
	// Time and Memory:

	static double cpuTime()
	{
		return (double)Stopwatch.GetTimestamp() / Stopwatch.Frequency;
	}

	static long memUsed()
	{
		return GC.GetTotalMemory(false);
	}

	[Conditional("DEBUG")]
	static public void assert(bool expr)
	{
		if (!expr)
			throw new Exception("assertion violated");
	}

	// Just like 'assert()' but expression will be evaluated in the release version as well.
	static void check(bool expr) { assert(expr); }

	// Redfine if you want output to go somewhere else:
	public static void reportf(string format, params object[]  args)
	{
		System.Console.Write(format, args);
	}
	public static void debug(string format, params object[]  args)
	{
		System.Console.WriteLine(format, args);
	}
	#endregion

    #region Stats, params
	public class SolverStats {
		public long   starts, decisions, propagations, conflicts;
		public long   clauses_literals, learnts_literals, max_literals, tot_literals;
	}


	public class SearchParams {
		public double  var_decay, clause_decay, random_var_freq;    // (reasonable values are: 0.95, 0.999, 0.02)    
		public SearchParams() : this(1,1,0) { }
		public SearchParams(SearchParams other) : this(other.var_decay, other.clause_decay, other.random_var_freq) { }
		public SearchParams(double v, double c, double r) { var_decay = v; clause_decay = c; 
		random_var_freq = r; }
	}
	#endregion

	#region VarOrder
	public class VarOrder {
		readonly protected vec<lbool>    assigns;     // var.val. Pointer to external assignment table.
		readonly protected vec<double>  activity;    // var.act. Pointer to external activity table.
		protected Heap   heap;
		double              random_seed; // For the internal random number generator

		public VarOrder(vec<lbool> ass, vec<double> act)
		{
			assigns = ass;
			activity = act;
			heap = new Heap(var_lt);
			random_seed = 91648253;
		}

		bool var_lt (Var x, Var y) { return activity[x] > activity[y]; }

		public virtual void newVar()
		{
			heap.setBounds(assigns.size());
			heap.insert(assigns.size()-1);
		}
		
		// Called when variable increased in activity.
		public virtual void update(Var x)
		{
			if (heap.inHeap(x))
				heap.increase(x);
		}

		// Called when variable is unassigned and may be selected again.
		public virtual void undo(Var x)
		{
			if (!heap.inHeap(x))
				heap.insert(x);
		}

		public Lit select()
		{
			return select(0.0);
		}

		// Selects a new, unassigned variable (or 'var_Undef' if none exists).
		public virtual Lit select(double random_var_freq)
		{
			// Random decision:
			if (drand(ref random_seed) < random_var_freq && !heap.empty()){
				Var next = irand(ref random_seed,assigns.size());
				if (isUndef(assigns[next]))
					return ~new Lit(next);
			}

			// Activity based decision:
			while (!heap.empty()){
				Var next = heap.getmin();
				if (isUndef(assigns[next]))
					return ~new Lit(next);	
			}

			return lit_Undef;
		}
	}	
	#endregion

	#region Solver state
	bool                ok;               // If FALSE, the constraints are already unsatisfiable. No part of the solver state may be used!
	protected vec<Clause>        clauses;          // List of problem clauses.
	protected vec<Clause>        learnts;          // List of learnt clauses.
	double              cla_inc;          // Amount to bump next clause with.
	double              cla_decay;        // INVERSE decay factor for clause activity: stores 1/decay.

	public vec<double>         activity;         // A heuristic measurement of the activity of a variable.
	double              var_inc;          // Amount to bump next variable with.
	double              var_decay;        // INVERSE decay factor for variable activity: stores 1/decay. Use negative value for static variable order.
	VarOrder            order;            // Keeps track of the decision variable order.

	vec<vec<Clause>>  watches;          // 'watches[lit]' is a list of constraints watching 'lit' (will go there if literal becomes true).
	public vec<lbool>           assigns;          // The current assignments.
	public vec<Lit>            trail;            // Assignment stack; stores all assigments made in the order they were made.
	protected vec<int>            trail_lim;        // Separator indices for different decision levels in 'trail'.
	protected vec<Clause>        reason;           // 'reason[var]' is the clause that implied the variables current value, or 'null' if none.
	protected vec<int>            level;            // 'level[var]' is the decision level at which assignment was made.
	vec<int>            trail_pos;        // 'trail_pos[var]' is the variable's position in 'trail[]'. This supersedes 'level[]' in some sense, and 'level[]' will probably be removed in future releases.
	int                 root_level;       // Level of first proper decision.
	int                 qhead;            // Head of queue (as index into the trail -- no more explicit propagation queue in MiniSat).
	int                 simpDB_assigns;   // Number of top-level assignments since last execution of 'simplifyDB()'.
	long               simpDB_props;     // Remaining number of propagations that must be made before next execution of 'simplifyDB()'.

	// Temporaries (to reduce allocation overhead). Each variable is prefixed by the method in which is used:
	//
	vec<lbool>           analyze_seen;
	vec<Lit>            analyze_stack;
	vec<Lit>            analyze_toclear;
	vec<Lit>            addUnit_tmp;
	vec<Lit>            addBinary_tmp;
	vec<Lit>            addTernary_tmp;
	#endregion

	#region Main internal methods:
	void        analyzeFinal     (Clause confl) { analyzeFinal(confl, false); }
	bool        enqueue          (Lit fact) { return enqueue(fact, null); }

	// Activity:
	//
	void     varBumpActivity(Lit p) {
		if (var_decay < 0) return;     // (negative decay means static variable order -- don't bump)
		if ( (activity[var(p)] += var_inc) > 1e100 ) varRescaleActivity();
		order.update(var(p)); }
	void     varDecayActivity  () { if (var_decay >= 0) var_inc *= var_decay; }
	void     claDecayActivity  () { cla_inc *= cla_decay; }

	// Operations on clauses:
	//
	void     newClause(vec<Lit> ps) { newClause(ps, false); }
	void     claBumpActivity (Clause c) { if ( (c.activity += (float)cla_inc) > 1e20 ) claRescaleActivity(); }
	protected void     remove          (Clause c) { remove(c, false); }
	protected bool     locked          (Clause c){ return c == reason[var(c[0])]; }

	int      decisionLevel() { return trail_lim.size(); }
	#endregion

	#region Public interface
	public Solver() {
				clauses = new vec<Clause>();
				learnts = new vec<Clause>();
				activity = new vec<double>();
				watches = new vec<vec<Clause>>();
				assigns = new vec<lbool>();
				trail_pos = new vec<int>();
				trail = new vec<Lit>();
				trail_lim = new vec<int>();
				reason = new vec<Clause>();
				level = new vec<int>();
				analyze_seen = new vec<lbool>();
				analyze_stack = new vec<Lit>();
				analyze_toclear = new vec<Lit>();
				addUnit_tmp = new vec<Lit>();
				addBinary_tmp = new vec<Lit>();
				addTernary_tmp = new vec<Lit>();
				model = new vec<lbool>();
				conflict = new vec<Lit>();
				addUnit_tmp .growTo(2);
				addBinary_tmp .growTo(2);
				addTernary_tmp.growTo(3);

				stats = new SolverStats();

			  ok               = true;
			  cla_inc          = 1;
			  cla_decay        = 1;
			  var_inc          = 1;
			  var_decay        = 1;
			  order            = createOrder();
			  qhead            = 0;
			  simpDB_assigns   = 0;
			  simpDB_props     = 0;
			  default_parms   = new SearchParams(0.95, 0.999, 0.02);
			  expensive_ccmin  = true;
			  verbosity        = 0;
			  progress_estimate= 0;
			 
				vec<Lit> dummy = new vec<Lit>(2,lit_Undef);
				dummy.pop();

			 }
	
	protected virtual VarOrder createOrder()
	{
		return new VarOrder(assigns, activity);
	}

   ~Solver() {
	   for (int i = 0; i < learnts.size(); i++) remove(learnts[i], true);
	   for (int i = 0; i < clauses.size(); i++) if (clauses[i] != null) remove(clauses[i], true); 

	}

	// Helpers: (semi-internal)
	//
	public lbool   value(Var x) { return assigns[x]; }
	public lbool   value(Lit p) { return sign(p) ? ~assigns[var(p)] : assigns[var(p)]; }

	public int     nAssigns() { return trail.size(); }
	public int     nClauses() { return clauses.size(); }   // (minor difference from MiniSat without the GClause trick: learnt binary clauses will be counted as original clauses)
	public int     nLearnts() { return learnts.size(); }

	// Statistics: (read-only member variable)
	//
	public SolverStats     stats;

	// Mode of operation:
	//
	public SearchParams    default_parms;     // Restart frequency etc.
	public bool            expensive_ccmin;    // Controls conflict clause minimization. TRUE by default.
	public int             verbosity;          // Verbosity level. 0=silent, 1=some progress report, 2=everything

	// Problem specification:
	//
	// public Var     newVar    ();
	public int     nVars     ()                    { return assigns.size(); }
	public void    addUnit   (Lit p)               { addUnit_tmp [0] = p; addClause (addUnit_tmp); }
	public void    addBinary (Lit p, Lit q)        { addBinary_tmp [0] = p; addBinary_tmp [1] = q; addClause(addBinary_tmp); }
	public void    addTernary(Lit p, Lit q, Lit r) { addTernary_tmp[0] = p; addTernary_tmp[1] = q; addTernary_tmp[2] = r; addClause(addTernary_tmp); }
	public void    addClause (vec<Lit> ps)  { newClause(ps); }  // (used to be a difference between internal and external method...)

	// Solving:
	//
	public bool    okay() { return ok; }       // FALSE means solver is in an conflicting state (must never be used again!)
	//public void    simplifyDB();
	//public bool    solve(vec<Lit> assumps);
	public bool    solve() { vec<Lit> tmp = new vec<Lit>(); return solve(tmp); }

	public double      progress_estimate;  // Set by 'search()'.
	public vec<lbool>  model;              // If problem is satisfiable, this vector contains the model (if any).
	public vec<Lit>    conflict;           // If problem is unsatisfiable (possibly under assumptions), this vector represent the conflict clause expressed in the assumptions.
	#endregion

	#region Operations on clauses:
	/*_________________________________________________________________________________________________
	|
	|  newClause : (ps : const vec<Lit>&) (learnt : bool)  .  [void]
	|  
	|  Description:
	|    Allocate and add a new clause to the SAT solvers clause database. If a conflict is detected,
	|    the 'ok' flag is cleared and the solver is in an unusable state (must be disposed).
	|  
	|  Input:
	|    ps     - The new clause as a vector of literals.
	|    learnt - Is the clause a learnt clause? For learnt clauses, 'ps[0]' is assumed to be the
	|             asserting literal. An appropriate 'enqueue()' operation will be performed on this
	|             literal. One of the watches will always be on this literal, the other will be set to
	|             the literal with the highest decision level.
	|  
	|  Effect:
	|    Activity heuristics are updated.
	|________________________________________________________________________________________________@*/
	vec<Lit> BasicClauseSimplification(vec<Lit> ps, bool copy)
	{
		vec<Lit> qs;
		if (copy) {
			qs = new vec<Lit>();
			ps.copyTo(qs);             // Make a copy of the input vector.
		} else {
			qs = ps;
		}

		Dictionary<Var, Lit> dict = new Dictionary<Var,Lit>(ps.size());
		int ptr = 0;

		for (int i = 0; i < qs.size(); i++) {
			Lit l = qs[i];
			Lit other;
			Var v = var(l);
			if (dict.TryGetValue(v, out other)) {
				if (other == l) {} // already seen it
				else return null; // other = ~l, so always satisfied
			} else {
				dict[v] = l;
				qs[ptr++] = l;
			}
		}
		qs.shrinkTo(ptr);

		return qs;
	}

	void reorderByLevel(vec<Lit> ps)
	{
		int max = int.MinValue;
		int max_at = -1;
		int max2 = int.MinValue;
		int max2_at = -1;
		for (int i = 0; i < ps.size(); ++i)
		{
			int lev = level[var(ps[i])];
			if (lev == -1) lev = int.MaxValue;
			else if (value(ps[i]) == lbool.True) lev = int.MaxValue;

			if (lev >= max) {
				max2_at = max_at;
				max2 = max;
				max = lev;
				max_at = i;
			} else if (lev > max2) {
				max2 = lev;
				max2_at = i;
			}
		}

		if (max_at == 0)
			ps.swap(1, max2_at);
		else if (max_at == 1)
			ps.swap(0, max2_at);
		else if (max2_at == 0)
			ps.swap(1, max_at);
		else if (max2_at == 1)
			ps.swap(0, max_at);
		else {
			ps.swap(0, max_at);
			ps.swap(1, max2_at);
		}
	}

	protected void newClause(vec<Lit> ps_, bool learnt)
	{
		newClause(ps_, learnt, false, true);
	}

	protected void newClause(vec<Lit> ps_, bool learnt, bool theoryClause, bool copy)
	{
		if (!ok) return;

		//foreach (Lit p in ps_) { Console.Write (" {0} ", p); } Console.WriteLine (" END");
		vec<Lit>    ps;

		assert(!(learnt && theoryClause));

		if (!learnt){
			assert(theoryClause || decisionLevel() == 0);

			vec<Lit> qs = BasicClauseSimplification(ps_, copy);

			if (qs == null) return;

			// Check if clause is satisfied:
			for (int i = 0; i < qs.size(); i++){
				if (level[var(qs[i])] == 0 && value(qs[i]) == l_True)
					return; }

			// Remove false literals:
			{
				int i, j; 
				for (i = j = 0; i < qs.size(); i++)
					if (level[var(qs[i])] != 0 || value(qs[i]) != l_False)
						qs[j++] = qs[i];
				qs.shrink(i - j);
			}

			ps = qs;
		} else
			ps = ps_;

		// 'ps' is now the (possibly) reduced vector of literals.

		if (ps.size() == 0){
			ok = false;

		}else if (ps.size() == 1){
			// NOTE: If enqueue takes place at root level, the assignment will be lost in incremental use (it doesn't seem to hurt much though).
			//if (!enqueue(ps[0], GClause_new(Clause_new(learnt, ps))))
			if (theoryClause) {
				levelToBacktrack = 0;
				cancelUntil(0);
			}

			Clause c = Clause_new(learnt || theoryClause, ps);
			NewClauseCallback(c);

			if (!enqueue(ps[0]))
				ok = false;
		}else{
			if (theoryClause)
				reorderByLevel(ps);

			// Allocate clause:
			Clause c   = Clause_new(learnt || theoryClause, ps);

			if (!learnt && !theoryClause) {
				// Store clause:
				clauses.push(c);
				stats.clauses_literals += c.size();
			} else {
				if (learnt) {
					// Put the second watch on the literal with highest decision level:
					int     max_i = 1;
					int     max   = level[var(ps[1])];
					for (int i = 2; i < ps.size(); i++)
						if (level[var(ps[i])] > max) {
							max   = level[var(ps[i])];
							max_i = i;
						}
					c[1]     = ps[max_i];
					c[max_i] = ps[1];

					check(enqueue(c[0], c));
				} else {
					MoveBack(c[0], c[1]);
				}

				// Bumping:
				claBumpActivity(c);         // (newly learnt clauses should be considered active)
				learnts.push(c);
				stats.learnts_literals += c.size();
			}


			// Watch clause:
			watches[index(~c[0])].push(c);
			watches[index(~c[1])].push(c);
			NewClauseCallback(c);
		}
	}


	// Disposes a clauses and removes it from watcher lists. NOTE! Low-level; does NOT change the 'clauses' and 'learnts' vector.
	//
	void remove(Clause c, bool just_dealloc)
	{
		if (!just_dealloc){
			removeWatch(watches[index(~c[0])], c);
			removeWatch(watches[index(~c[1])], c);
		}

		if (c.learnt()) stats.learnts_literals -= c.size();
		else             stats.clauses_literals -= c.size();

		//xfree(c);
	}


	// Can assume everything has been propagated! (esp. the first two literals are != l_False, unless
	// the clause is binary and satisfied, in which case the first literal is true)
	// Returns True if clause is satisfied (will be removed), False otherwise.
	//
	bool simplify(Clause c)
	{
		assert(decisionLevel() == 0);
		for (int i = 0; i < c.size(); i++){
			if (value(c[i]) == l_True)
				return true;
		}
		return false;
	}
	#endregion

	#region Minor methods
	static bool removeWatch(vec<Clause> ws, Clause elem)    // Pre-condition: 'elem' must exists in 'ws' OR 'ws' must be empty.
	{
		if (ws.size() == 0) return false;     // (skip lists that are already cleared)
		int j = 0;
		for (; ws[j] != elem  ; j++) assert(j < ws.size() - 1);
		for (; j < ws.size()-1; j++) ws[j] = ws[j+1];
		ws.pop();
		return true;
	}
	// Creates a new SAT variable in the solver. If 'decision_var' is cleared, variable will not be
	// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
	//
	public Var newVar() {
		int     index;
		index = nVars();
		watches     .push(new vec<Clause>());          // (list for positive literal)
		watches     .push(new vec<Clause>());          // (list for negative literal)
		reason      .push(null);
		assigns     .push(lbool.Undef0);
		level       .push(-1);
		trail_pos	.push(-1); 
		activity    .push(0);
		order       .newVar();
		analyze_seen.push(0);
		return index; }


	// Returns FALSE if immediate conflict.
	bool assume(Lit p) {
		trail_lim.push(trail.size());
		return enqueue(p); }


	// Revert to the state at given level.
	protected void cancelUntil(int level) {
		CancelUntilCallback(level);
		if (decisionLevel() > level){
			for (int c = trail.size()-1; c >= trail_lim[level]; c--){
				Var     x  = var(trail[c]);
				assigns[x] = lbool.Undef0;
				reason [x] = null;
				order.undo(x); }
			trail.shrink(trail.size() - trail_lim[level]);
			trail_lim.shrink(trail_lim.size() - level);
			qhead = trail.size(); } }
	#endregion

	#region Major methods:
	/*_________________________________________________________________________________________________
	|
	|  analyze : (confl : Clause*) (out_learnt : vec<Lit>&) (out_btlevel : int&)  .  [void]
	|  
	|  Description:
	|    Analyze conflict and produce a reason clause.
	|  
	|    Pre-conditions:
	|      * 'out_learnt' is assumed to be cleared.
	|      * Current decision level must be greater than root level.
	|  
	|    Post-conditions:
	|      * 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.
	|  
	|  Effect:
	|    Will undo part of the trail, upto but not beyond the assumption of the current decision level.
	|________________________________________________________________________________________________@*/
	void analyze(Clause confl, vec<Lit> out_learnt, out int out_btlevel)
	{
		vec<lbool>     seen  = analyze_seen;
		int            pathC = 0;
		Lit            p     = lit_Undef;

		AdditionalConflictAnalisis(confl.GetData(), confl);

		// Generate conflict clause:
		//
		out_learnt.push(new Lit());      // (leave room for the asserting literal)
		out_btlevel = 0;
		int index = trail.size()-1;
		//debug("start analyze");
		do{
		/*
		    debug("    loop analyze {0} {1} {2}\n", confl, p, p==lit_Undef ? -1 : level[var(p)]);
			if (confl == null)
			{
			  for (int i = trail.size()-1; i >= 0; i--)
			    debug("   {0} {1} {2} {3}\n", trail[i], seen[var(trail[i])], 
						level[var(trail[i])], reason[var(trail[i])]);
			} */
			assert(confl != null);          // (otherwise should be UIP)

			Clause c = confl;

			if (c.learnt())
				claBumpActivity(c);

			for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
				Lit q = c[j];
				if (seen[var(q)] == 0 && level[var(q)] > 0){
					varBumpActivity(q);
					seen[var(q)] = lbool.True;
					if (level[var(q)] == decisionLevel())
						pathC++;
					else{
						out_learnt.push(q);
						out_btlevel = Math.Max(out_btlevel, level[var(q)]);
					}
				}
			}

			// Select next clause to look at:
			while (seen[var(trail[index--])] == 0);
			p     = trail[index+1];
			confl = reason[var(p)];
			seen[var(p)] = 0;
			pathC--;

		}while (pathC > 0);
		out_learnt[0] = ~p;

		// Conflict clause minimization:
		{
		int     i, j;
		if (expensive_ccmin){
			// Simplify conflict clause (a lot):
			//
			uint    min_level = 0;
			for (i = 1; i < out_learnt.size(); i++)
				min_level |= (uint)(1 << (level[var(out_learnt[i])] & 31));         // (maintain an abstraction of levels involved in conflict)

			analyze_toclear.clear();
			for (i = j = 1; i < out_learnt.size(); i++)
				if (reason[var(out_learnt[i])] == null || !analyze_removable(out_learnt[i], min_level))
					out_learnt[j++] = out_learnt[i];
		}else{
			// Simplify conflict clause (a little):
			//
			analyze_toclear.clear();
			for (i = j = 1; i < out_learnt.size(); i++){
				Clause r = reason[var(out_learnt[i])];
				if (r == null)
					out_learnt[j++] = out_learnt[i];
				else{
					Clause c = r;
					for (int k = 1; k < c.size(); k++)
						if (seen[var(c[k])]==0 && level[var(c[k])] != 0){
							out_learnt[j++] = out_learnt[i];
							goto Keep; }
					analyze_toclear.push(out_learnt[i]);
					Keep: ;
				}
			}
		}


		// Clean up:
		//
		{
		int jj;
		for (jj = 0; jj < out_learnt.size()     ; jj++) seen[var(out_learnt     [jj])] = 0;
		for (jj = 0; jj < analyze_toclear.size(); jj++) seen[var(analyze_toclear[jj])] = 0;    // ('seen[]' is now cleared)
		}


		stats.max_literals += out_learnt.size();
		out_learnt.shrink(i - j);
		stats.tot_literals += out_learnt.size();

		}
	}


	// Check if 'p' can be removed. 'min_level' is used to abort early if visiting literals at a level that cannot be removed.
	//
	bool analyze_removable(Lit p_, uint min_level)
	{
		assert(reason[var(p_)] != null);
		analyze_stack.clear(); analyze_stack.push(p_);
		int top = analyze_toclear.size();
		while (analyze_stack.size() > 0){
			assert(reason[var(analyze_stack.last())] != null);
			Clause c = reason[var(analyze_stack.last())]; 
			analyze_stack.pop();
			for (int i = 1; i < c.size(); i++){
				Lit p = c[i];
				if (analyze_seen[var(p)]==0 && level[var(p)] != 0){
					if (reason[var(p)] != null && ((1 << (level[var(p)] & 31)) & min_level) != 0){
						analyze_seen[var(p)] = lbool.True;
						analyze_stack.push(p);
						analyze_toclear.push(p);
					}else{
						for (int j = top; j < analyze_toclear.size(); j++)
							analyze_seen[var(analyze_toclear[j])] = 0;
						analyze_toclear.shrink(analyze_toclear.size() - top);
						return false;
					}
				}
			}
		}

		analyze_toclear.push(p_);

		return true;
	}


	/*_________________________________________________________________________________________________
	|
	|  analyzeFinal : (confl : Clause*) (skip_first : bool)  .  [void]
	|  
	|  Description:
	|    Specialized analysis procedure to express the final conflict in terms of assumptions.
	|    'root_level' is allowed to point beyond end of trace (useful if called after conflict while
	|    making assumptions). If 'skip_first' is TRUE, the first literal of 'confl' is  ignored (needed
	|    if conflict arose before search even started).
	|________________________________________________________________________________________________@*/
	void analyzeFinal(Clause confl, bool skip_first)
	{
		// -- NOTE! This code is relatively untested. Please report bugs!
		conflict.clear();
		if (root_level == 0) return;

		vec<lbool>     seen  = analyze_seen;
		for (int i = skip_first ? 1 : 0; i < confl.size(); i++){
			Var     x = var(confl[i]);
			if (level[x] > 0)
				seen[x] = lbool.True;
		}

		int     start = (root_level >= trail_lim.size()) ? trail.size()-1 : trail_lim[root_level];
		for (int i = start; i >= trail_lim[0]; i--){
			Var     x = var(trail[i]);
			if (seen[x]!=0){
				Clause r = reason[x];
				if (r == null){
					assert(level[x] > 0);
					conflict.push(~trail[i]);
				}else{
					Clause c = r;
					for (int j = 1; j < c.size(); j++)
						if (level[var(c[j])] > 0)
							seen[var(c[j])] = lbool.True;
				}
				seen[x] = lbool.Undef0;
			}
		}
	}


	/*_________________________________________________________________________________________________
	|
	|  enqueue : (p : Lit) (from : Clause*)  .  [bool]
	|  
	|  Description:
	|    Puts a new fact on the propagation queue as well as immediately updating the variable's value.
	|    Should a conflict arise, FALSE is returned.
	|  
	|  Input:
	|    p    - The fact to enqueue
	|    from - [Optional] Fact propagated from this (currently) unit clause. Stored in 'reason[]'.
	|           Default value is null (no reason).
	|  
	|  Output:
	|    TRUE if fact was enqueued without conflict, FALSE otherwise.
	|________________________________________________________________________________________________@*/
	bool enqueue(Lit p, Clause from)
	{
		if (!isUndef(value(p))) {
			return value(p) != l_False;
		}else{
			Var x = var(p);
			assigns[x] = toLbool(!sign(p));
			level  [x] = decisionLevel();
			trail_pos[x] = trail.size();
			reason [x] = from;
			trail.push(p);
			return true;
		}
	}


	/*_________________________________________________________________________________________________
	|
	|  propagate : [void]  .  [Clause*]
	|  
	|  Description:
	|    Propagates all enqueued facts. If a conflict arises, the conflicting clause is returned,
	|    otherwise null. NOTE! This method has been optimized for speed rather than readability.
	|  
	|    Post-conditions:
	|      * the propagation queue is empty, even if there was a conflict.
	|________________________________________________________________________________________________@*/
	Clause propagate()
	{
		Clause confl = null;
		while (qhead < trail.size()){
			stats.propagations++;
			simpDB_props--;

			Lit            p   = trail[qhead++];     // 'p' is enqueued fact to propagate.
			vec<Clause>    ws  = watches[index(p)];
			//GClause*       i,* j, *end;
			int i, j, end;

			for (i = j = 0, end = i + ws.size();  i != end;){
				Clause c = ws[i++];
				// Make sure the false literal is data[1]:
				Lit false_lit = ~p;
				if (c[0] == false_lit)
					{ c[0] = c[1]; c[1] = false_lit; }

				assert(c[1] == false_lit);

				// If 0th watch is true, then clause is already satisfied.
				Lit   first = c[0];
				lbool val   = value(first);
				if (val == l_True){
					ws[j++] = c;
				}else{
					// Look for new watch:
					for (int k = 2; k < c.size(); k++)
						if (value(c[k]) != l_False){
							c[1] = c[k]; c[k] = false_lit;
							watches[index(~c[1])].push(c);
							goto FoundWatch; }

					// Did not find watch -- clause is unit under assignment:
					ws[j++] = c;
					if (!enqueue(first, c)){
						if (decisionLevel() == 0)
							ok = false;
						confl = c;
						qhead = trail.size();
						// Copy the remaining watches:
						while (i < end)
							ws[j++] = ws[i++];
					}
				  FoundWatch:;
				}
			}
			ws.shrink(i - j);
		}

		return confl;
	}


	/*_________________________________________________________________________________________________
	|
	|  reduceDB : ()  .  [void]
	|  
	|  Description:
	|    Remove half of the learnt clauses, minus the clauses locked by the current assignment. Locked
	|    clauses are clauses that are reason to some assignment. Binary clauses are never removed.
	|________________________________________________________________________________________________@*/
	class reduceDB_lt : IComparer<Clause>
	{ 
	  public int Compare(Clause x, Clause y) { 
		if (x.size() > 2 && (y.size() == 2 || x.activity < y.activity)) 
		return -1; else return 1; 
	  } 
	}

	void reduceDB()
	{
		int     i, j;
		double  extra_lim = cla_inc / learnts.size();    // Remove any clause below this activity

		learnts.sort(new reduceDB_lt());
		for (i = j = 0; i < learnts.size() / 2; i++){
			if (learnts[i].size() > 2 && !locked(learnts[i]))
				remove(learnts[i]);
			else
				learnts[j++] = learnts[i];
		}
		for (; i < learnts.size(); i++){
			if (learnts[i].size() > 2 && !locked(learnts[i]) && learnts[i].activity < extra_lim)
				remove(learnts[i]);
			else
				learnts[j++] = learnts[i];
		}
		learnts.shrink(i - j);
	}


	/*_________________________________________________________________________________________________
	|
	|  simplifyDB : [void]  .  [bool]
	|  
	|  Description:
	|    Simplify the clause database according to the current top-level assigment. Currently, the only
	|    thing done here is the removal of satisfied clauses, but more things can be put here.
	|________________________________________________________________________________________________@*/
	void simplifyDB()
	{
		if (!ok) return;    // GUARD (public method)
		assert(decisionLevel() == 0);

		if (propagate() != null){
			ok = false;
			return; }

		if (nAssigns() == simpDB_assigns || simpDB_props > 0)   // (nothing has changed or preformed a simplification too recently)
			return;

		// Clear watcher lists:
		for (int i = simpDB_assigns; i < nAssigns(); i++){
			Lit           p  = trail[i];
			watches[index( p)].clear();
			watches[index(~p)].clear();
		}

		// Remove satisfied clauses:
		for (int type = 0; type < 2; type++){
			vec<Clause> cs = type != 0 ? learnts : clauses;
			int     j  = 0;
			for (int i = 0; i < cs.size(); i++){
				if (!locked(cs[i]) && simplify(cs[i])) {
					remove(cs[i]);
					RemoveClauseCallback(cs[i]);
				} else
					cs[j++] = cs[i];
			}
			cs.shrink(cs.size()-j);
		}

		simpDB_assigns = nAssigns();
		simpDB_props   = stats.clauses_literals + stats.learnts_literals;   // (shouldn't depend on 'stats' really, but it will do for now)
	}


	/*_________________________________________________________________________________________________
	|
	|  search : (nof_conflicts : int) (nof_learnts : int) (parms : const SearchParams&)  .  [lbool]
	|  
	|  Description:
	|    Search for a model the specified number of conflicts, keeping the number of learnt clauses
	|    below the provided limit. NOTE! Use negative value for 'nof_conflicts' or 'nof_learnts' to
	|    indicate infinity.
	|  
	|  Output:
	|    'l_True' if a partial assigment that is consistent with respect to the clauseset is found. If
	|    all variables are decision variables, this means that the clause set is satisfiable. 'l_False'
	|    if the clause set is unsatisfiable. 'l_Undef' if the bound on number of conflicts is reached.
	|________________________________________________________________________________________________@*/
	lbool search(int nof_conflicts, int nof_learnts, SearchParams parms)
	{
		if (!ok) return l_False;    // GUARD (public method)
		assert(root_level == decisionLevel());

		stats.starts++;
		int     conflictC = 0;
		var_decay = 1 / parms.var_decay;
		cla_decay = 1 / parms.clause_decay;
		model.clear();

		for (;;){
			Clause confl = propagate();
			if (confl != null){
				// CONFLICT

				stats.conflicts++; conflictC++;
				vec<Lit>    learnt_clause = new vec<Lit>();
				int         backtrack_level;
				if (decisionLevel() == root_level){
					// Contradiction found:
					analyzeFinal(confl);
					return l_False; }
				analyze(confl, learnt_clause, out backtrack_level);
				cancelUntil(Math.Max(backtrack_level, root_level));
				newClause(learnt_clause, true);
				if (learnt_clause.size() == 1) level[var(learnt_clause[0])] = 0;    // (this is ugly (but needed for 'analyzeFinal()') -- in future versions, we will backtrack past the 'root_level' and redo the assumptions)
				varDecayActivity();
				claDecayActivity();

			}else{
				// NO CONFLICT

				if (nof_conflicts >= 0 && conflictC >= nof_conflicts){
					// Reached bound on number of conflicts:
					progress_estimate = progressEstimate();
					cancelUntil(root_level);
					return lbool.Undef0; }

				if (decisionLevel() == 0)
					// Simplify the set of problem clauses:
					{ simplifyDB(); if (!ok) return l_False; }

				if (nof_learnts >= 0 && learnts.size()-nAssigns() >= nof_learnts)
					// Reduce the set of learnt clauses:
					reduceDB();

				// New variable decision:
				stats.decisions++;
				Lit next = order.select(parms.random_var_freq);

				if (next == lit_Undef){
				    if (ModelFound())
						continue;
					// Model found:
					model.growTo(nVars());
					for (int i = 0; i < nVars(); i++) model[i] = value(i);
					cancelUntil(root_level);
					return l_True;
				}

				check(assume(next));
			}
		}
	}


	// Divide all variable activities by 1e100.
	//
	void varRescaleActivity()
	{
		for (int i = 0; i < nVars(); i++)
			activity[i] *= 1e-100;
		var_inc *= 1e-100;
	}


	// Divide all constraint activities by 1e20.
	//
	void claRescaleActivity()
	{
		for (int i = 0; i < learnts.size(); i++)
			learnts[i].activity *= 1e-20f;
		cla_inc *= 1e-20;
	}


	/*_________________________________________________________________________________________________
	|
	|  solve : (assumps : const vec<Lit>&)  .  [bool]
	|  
	|  Description:
	|    Top-level solve. If using assumptions (non-empty 'assumps' vector), you must call
	|    'simplifyDB()' first to see that no top-level conflict is present (which would put the solver
	|    in an undefined state).
	|
	|  Input:
	|    A list of assumptions (unit clauses coded as literals). Pre-condition: The assumptions must
	|    not contain both 'x' and '~x' for any variable 'x'.
	|________________________________________________________________________________________________@*/
	bool solve(vec<Lit> assumps)
	{
		simplifyDB();
		if (!ok) return false;

		SearchParams    parms = new SearchParams(default_parms);
		double  nof_conflicts = 100;
		double  nof_learnts   = nClauses() / 3;
		lbool   status        = lbool.Undef0;

		// Perform assumptions:
		root_level = assumps.size();
		for (int i = 0; i < assumps.size(); i++){
			Lit p = assumps[i];
			assert(var(p) < nVars());
			if (!assume(p)){
				Clause r = reason[var(p)];
				if (r != null){
					analyzeFinal(r, true);
					conflict.push(~p);
				}else {
					conflict.clear();
					conflict.push(~p);
				}
				cancelUntil(0);
				return false; 
			}

			{
				Clause confl = propagate();
				if (confl != null){
					analyzeFinal(confl);
					assert(conflict.size() > 0);
					cancelUntil(0);
					return false; 
				}
			}
		}
		assert(root_level == decisionLevel());

		// Search:
		if (verbosity >= 1){
			reportf("==================================[MINISAT]===================================\n");
			reportf("| Conflicts |     ORIGINAL     |              LEARNT              | Progress |\n");
			reportf("|           | Clauses Literals |   Limit Clauses Literals  Lit/Cl |          |\n");
			reportf("==============================================================================\n");
		}

		while (isUndef(status)){
			if (verbosity >= 1)
				reportf("| {0,9} | {1,7} {2,8} | {3,7} {4,7} {5,8} {6,7:0.0} |{7,6:0.000} %% |\n",
					(int)stats.conflicts, nClauses(), (int)stats.clauses_literals, 
					(int)nof_learnts, nLearnts(), (int)stats.learnts_literals, 
					(double)stats.learnts_literals/nLearnts(), progress_estimate*100);
			status = search((int)nof_conflicts, (int)nof_learnts, parms);
			nof_conflicts *= 1.5;
			nof_learnts   *= 1.1;
		}
		if (verbosity >= 1)
			reportf("==============================================================================\n");

		cancelUntil(0);
		return status == l_True;
	}
	#endregion

	#region Stats
	double start_time = cpuTime();

	// Return search-space coverage. Not extremely reliable.
	//
	double progressEstimate()
	{
		double  progress = 0;
		double  F = 1.0 / nVars();
		for (int i = 0; i < nVars(); i++)
			if (!isUndef(value(i)))
				progress += Math.Pow(F, level[i]);
		return progress / nVars();
	}


	public void printStats()
	{
		double  cpu_time = cpuTime()-start_time;
		long   mem_used = memUsed();
		reportf("restarts              : {0,12}\n", stats.starts);
		reportf("conflicts             : {0,12}   ({1:0.0} /sec)\n", stats.conflicts   , stats.conflicts   /cpu_time);
		reportf("decisions             : {0,12}   ({1:0.0} /sec)\n", stats.decisions   , stats.decisions   /cpu_time);
		reportf("propagations          : {0,12}   ({1:0.0} /sec)\n", stats.propagations, stats.propagations/cpu_time);
		reportf("conflict literals     : {0,12}   ({1:0.00} %% deleted)\n", stats.tot_literals, (stats.max_literals - stats.tot_literals)*100 / (double)stats.max_literals);
		if (mem_used != 0) reportf("Memory used           : {0:0.00} MB\n", mem_used / 1048576.0);
		reportf("CPU time              : {0:0.000} s\n", cpu_time);
	}
	#endregion


	#region DPLL(T) stuff -- don't use yet
	// cancelUntil was called with level
	protected virtual void CancelUntilCallback(int level)
	{
	}

	int levelToBacktrack;

	bool ModelFound()
	{
		levelToBacktrack = int.MaxValue;

		bool res = ModelFoundCallback();

		if (levelToBacktrack != int.MaxValue) {
			cancelUntil(levelToBacktrack);
			qhead = trail_lim.size() == 0 ? 0 : trail_lim.last();
		}

		if (!ok)
			return false;
		
		return res;
	}

	void MoveBack(Lit l1, Lit l2)
	{
		
		int lev1 = level[var(l1)];
		int lev2 = level[var(l2)];
		if (lev1 == -1) lev1 = int.MaxValue;
		if (lev2 == -1) lev2 = int.MaxValue;
		if (lev1 < levelToBacktrack || lev2 < levelToBacktrack) {
			if (value(l1) == lbool.True) {
				if (value(l2) == lbool.True)
					{}
				else if (lev1 <= lev2 || levelToBacktrack <= lev2)
					{}
				else
					levelToBacktrack = lev2;
			} else {
				if (value(l2) == lbool.True) {
					if (lev2 <= lev1 || levelToBacktrack <= lev1)
						{}
					else
						levelToBacktrack = lev1;
				} else
					levelToBacktrack = Math.Min (lev1, lev2);
			}
					
		}

		//debug("level: {0} {1}", levelToBacktrack, l);
	}

	// this is expected to return true if it adds some new conflict clauses
	protected virtual bool ModelFoundCallback()
	{
		return false;
	}

	public bool SearchNoRestarts()
	{
		SearchParams parms = new SearchParams(default_parms);
		lbool status = lbool.Undef0;

		simplifyDB();
		if (!ok) return false;

		root_level = 0;
		assert(root_level == decisionLevel());

		while (isUndef(status)){
			status = search(-1, -1, parms);
		}

		cancelUntil(0);
		return ok && status == l_True;
	}

	protected virtual void NewClauseCallback(Clause c)
	{
	}

	protected virtual void RemoveClauseCallback(Clause c)
	{
	}

    protected virtual void AdditionalConflictAnalisis(Lit[] conflict, Clause c)
	{
	}
	#endregion

} // end class Solver

} // end namespace
