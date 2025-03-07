/******************************************************************************

Copyright (c) 2010, Cormac Flanagan (University of California, Santa Cruz)
                    and Stephen Freund (Williams College) 

All rights reserved.  

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

 * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

 * Neither the names of the University of California, Santa Cruz
      and Williams College nor the names of its contributors may be
      used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 ******************************************************************************/

package tools.hb;

import java.io.PrintWriter;
import java.util.Random;
import java.util.Vector;
import java.util.concurrent.SynchronousQueue;
import java.util.*;

import rr.annotations.Abbrev;
import rr.barrier.BarrierEvent;
import rr.barrier.BarrierListener;
import rr.barrier.BarrierMonitor;
import rr.error.ErrorMessage;
import rr.error.ErrorMessages;
import rr.event.AccessEvent;
import rr.event.AcquireEvent;
import rr.event.ArrayAccessEvent;
import rr.event.FieldAccessEvent;
import rr.event.InterruptEvent;
import rr.event.InterruptedEvent;
import rr.event.JoinEvent;
import rr.event.NewThreadEvent;
import rr.event.NotifyEvent;
import rr.event.ReleaseEvent;
import rr.event.StartEvent;
import rr.event.VolatileAccessEvent;
import rr.event.WaitEvent;
import rr.event.AccessEvent.Kind;
import rr.meta.ArrayAccessInfo;
import rr.meta.FieldInfo;
import rr.state.ShadowLock;
import rr.state.ShadowThread;
import rr.state.ShadowVar;
import rr.tool.Tool;
import tools.util.VectorClock;
import tools.util.VectorClockPair;
import acme.util.Assert;
import acme.util.Util;
import acme.util.decorations.Decoration;
import acme.util.decorations.DecorationFactory;
import acme.util.decorations.DefaultValue;
import acme.util.decorations.NullDefault;
import acme.util.option.CommandLine;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
/**
 * A simple VC-based HappensBefore Race Detector.
 *
 * This does not handle many special cases related to static initializers, etc.
 * and may report spurious warnings as a result.  The FastTrack implementations
 * do handles those items.
 */

@Abbrev("HB")
public final class HappensBeforeTool extends Tool implements BarrierListener<HBBarrierState> {

	/* Reporters for field/array errors */
	public final ErrorMessage<FieldInfo> errors = ErrorMessages.makeFieldErrorMessage("HappensBefore");
	public final ErrorMessage<ArrayAccessInfo> arrayErrors = ErrorMessages.makeArrayErrorMessage("HappensBefore");

	// my stuff
	//private PrintWriter mylog; 
	BufferedWriter mylog;
	private Integer AccessCounter = new Integer(0);
	private Integer LockCounter  = new Integer(0);
	private Integer SigCounter = new Integer(0);
	private Random idgen = new Random(); 
	private HashMap hm = new HashMap();
	private Integer locCount = new Integer(0);

	public HappensBeforeTool(String name, Tool next, CommandLine commandLine) {
		super(name, next, commandLine); 

		/* 
		 * Create a barrier monitor that will notify this tool when barrier ops happen.
		 * The get(k) method returns the initial state to associate with each barrier
		 * when the barrier is created.
		 */
		new BarrierMonitor<HBBarrierState>(this, new DefaultValue<Object,HBBarrierState>() {
			public HBBarrierState get(Object k) {
				return new HBBarrierState(ShadowLock.get(k));
			}
		});

		try {
			//mylog = new PrintWriter("my.log", "UTF-8");
            mylog = new BufferedWriter(new FileWriter("/mnt/my.log"), 32768);
		} catch (Exception e) {
			System.out.println("bad!");
		}

	}


	/* 
	 * Special methods that tells the instrumentor to create a field in ShadowThread
	 * named "hb" of type "VectorClock".  This is for performance only -- you could do the same
	 * thing using a decoration on ThreadStates.
	 */
	static VectorClock ts_get_cv_hb(ShadowThread ts) { Assert.panic("Bad");	return null; }
	static void ts_set_cv_hb(ShadowThread ts, VectorClock cv) { Assert.panic("Bad");  }

	private VectorClock get(ShadowThread td) {
		return ts_get_cv_hb(td);
	}

	/*
	 * Attach a VectorClock to each object used as a lock. 
	 */
	Decoration<ShadowLock,VectorClock> shadowLock = ShadowLock.decoratorFactory.make("HB:lock", DecorationFactory.Type.MULTIPLE,
			new DefaultValue<ShadowLock,VectorClock>() { public VectorClock get(ShadowLock ld) { return new VectorClock(1); }});

	private VectorClock get(ShadowLock td) {
		return shadowLock.get(td);
	}



	@Override
	public void create(NewThreadEvent e) {
		ShadowThread td = e.getThread();
		VectorClock cv = new VectorClock(td.getTid() + 1);
		ts_set_cv_hb(td, cv);
		super.create(e);

	}

	/*
	 * increment the time of the current thread.
	 */
	public void tick(ShadowThread currentThread) {
		get(currentThread).tick(currentThread.getTid());
	}

	@Override
	public void acquire(AcquireEvent ae) {
		final ShadowThread currentThread = ae.getThread();
		final ShadowLock shadowLock = ae.getLock();

	
		

		tick(currentThread);
		synchronized(shadowLock) {
			get(currentThread).max(get(shadowLock));
			synchronized (mylog) {
				++LockCounter;
				try {
					mylog.write((currentThread.getTid()+1) + ",LK," +  (shadowLock.hashCode()+1) + ",nil," + LockCounter +"\n");
					//mylog.write((currentThread.getTid()+1) + ",LK," +  (shadowLock.hashCode()+1) + ",nil," + LockCounter+"\n");
				} catch (Exception e) {
					System.out.println("bad!");
				}
				//mylog.flush();
			}
		}
		super.acquire(ae);
	}

	@Override
	public void release(ReleaseEvent re) {
		final ShadowThread currentThread = re.getThread();
		final ShadowLock shadowLock = re.getLock();

		

		synchronized(shadowLock) {
			get(shadowLock).copy(get(currentThread));
			synchronized (mylog) {
				//mylog.write("UNLOCK,T " +  currentThread.getTid() + ",VAR " +  shadowLock.hashCode());
					try {
						mylog.write((currentThread.getTid()+1) + ",UK," +  (shadowLock.hashCode()+1) + ",nil,0\n");
					//	mylog.write((currentThread.getTid()+1) + ",UK," +  (shadowLock.hashCode()+1) + ",nil,0\n");
					} catch (Exception e) {
						System.out.println("bad!");
					}
				//mylog.flush();
			}
		}
		tick(currentThread);
		super.release(re);
	}

	@Override
	public void volatileAccess(VolatileAccessEvent fae) {
		ShadowVar g = fae.getOriginalShadow();
		final ShadowThread currentThread = fae.getThread();

		if (g instanceof VectorClockPair) {
			final ShadowVar orig = fae.getOriginalShadow();
			final ShadowThread td = fae.getThread();
			VectorClockPair p = (VectorClockPair)g;

			final VectorClock cv = ts_get_cv_hb(td);
			Integer loc = new Integer(0);
			String s = fae.getAccessInfo().getLoc().toString();
			if (hm.containsKey(s)) {
				loc = (Integer)hm.get(s);
			} else {
				loc = ++locCount;
				hm.put(s, loc);
			}
			if (fae.isWrite()) {
				// synchronized (mylog) {
				// 	p.LastThread = currentThread.getTid()+1;
				// 	++AccessCounter;
				// 	try {				
				// 		mylog.write((currentThread.getTid()+1) + ",AWR," +  p.Identity + "," +  s + "," + AccessCounter+"\n");
				// 	} catch (Exception e) {
				// 		System.out.println("bad!");
				// 	}
				// //mylog.write("WRITE,T " +  tid + ",VAR " +  p.Identity + ", LOC " + fae.getAccessInfo().getLoc());
				// 	//mylog.flush();
				// }
				p.rd.max(get(currentThread));
				tick(td); 			
			} else {
				// synchronized (mylog) {
				// 	p.LastThread = currentThread.getTid()+1;
				// 	++AccessCounter;
				// 	try {
				// 		mylog.write((currentThread.getTid()+1) + ",ARD," +  p.Identity + "," +  s + "," + AccessCounter+"\n");
				// 	} catch (Exception e) {
				// 		System.out.println("bad!");
				// 	}
				// //mylog.write("WRITE,T " +  tid + ",VAR " +  p.Identity + ", LOC " + fae.getAccessInfo().getLoc());
				// 	//mylog.flush();
				// }
				synchronized(p.rd) {
					get(currentThread).max(p.rd);
				}
			}
		}
		super.volatileAccess(fae);

	}

	// FIXME: syncrhonize-me because we have destructive races on clock-vectors, and
	// we have no "right-mover" argument
	@Override
	public void access(AccessEvent fae) {

		ShadowVar g = fae.getOriginalShadow();
		final ShadowThread currentThread = fae.getThread();


		if (g instanceof VectorClockPair) {
			boolean passAlong = false;
			VectorClockPair p = (VectorClockPair)g;
			boolean isWrite = fae.isWrite();
			VectorClock cv = get(currentThread);
//			Util.log("p=" + p);
//			Util.log("t=" + cv);
			final int tid = currentThread.getTid();
			Object target = fae.getTarget();
			
			Integer loc = new Integer(0);
			String s = fae.getAccessInfo().getLoc().toString();
			if (hm.containsKey(s)) {
				loc = (Integer)hm.get(s);
			} else {
				loc = ++locCount;
				hm.put(s, loc);
			}

			if (p.LastThread == 0) {
				p.LastThread = tid; //init
			}

			if (!p.isShared && p.LastThread != tid) {
				p.isShared = true;
				synchronized(mylog) {
					try {
				
						mylog.write(p.lastmsg);

					} catch (Exception e) {
						System.out.println("bad!");
					}
				}
			}
			
			if (isWrite) {
				
				// check after prev read
				passAlong |= checkAfter(p.rd, "read", currentThread, "write", fae, true, p);
				// check after prev write
				passAlong |= checkAfter(p.wr, "write", currentThread, "write", fae, true, p);
				synchronized(p.wr) { 	
					p.wr.set(tid, cv.get(tid));
					synchronized (mylog) {
						++AccessCounter;
						try {
							if (p.isShared) {
								mylog.write((currentThread.getTid()+1) + ",WR," +  p.Identity + "," +  s + "," + AccessCounter+"\n");
							} else {
								p.lastmsg = (currentThread.getTid()+1) + ",WR," +  p.Identity + "," +  s + "," + AccessCounter+"\n";
							}
						} catch (Exception e) {
							System.out.println("bad!");
						}
					}
				}
 
			} else {
				
				// check after prev write
				passAlong |= checkAfter(p.wr, "write", currentThread, "read", fae, true, p);
				synchronized(p.rd) { 	
					p.rd.set(tid, cv.get(tid));
					synchronized (mylog) {
						++AccessCounter;
						try {
							if (p.isShared) {
								mylog.write((currentThread.getTid()+1) + ",RD," +  p.Identity + "," +  s + "," + AccessCounter+"\n");
							} else {
								p.lastmsg = (currentThread.getTid()+1) + ",RD," +  p.Identity + "," +  s + "," + AccessCounter+"\n";
							}
						} catch (Exception e) {
							System.out.println("bad!");
						}
					//mylog.write("READ,T " +  tid + ",VAR " +  p.Identity + ", LOC " + fae.getAccessInfo().getLoc());
						//mylog.flush();
					}
				}
//				p.rd.max(cv);
			}
			if (passAlong) {
				advance(fae);
			}
//			Util.log("p=" + p);
//			Util.log("t=" + cv);
		} else {
			super.access(fae);
		} 
	}
	
	// public static boolean readFastPath(final ShadowVar shadow, final ShadowThread st) {		
	// 	if (shadow instanceof VectorClockPair) {
	// 		VectorClockPair p = (VectorClockPair)shadow;

	// 		if (p.LastThread == st.getTid()) {
	// 			p.isShared = true;
	// 		}
	// 	}
	// 	return false;
	// }
	
	// public static boolean writeFastPath(final ShadowVar shadow, final ShadowThread st) {
	// 	if (shadow instanceof VectorClockPair) {
	// 		VectorClockPair p = (VectorClockPair)shadow;
	// 		if (p.LastThread == st.getTid()) {
	// 			return true;
	// 		}
	// 	}
	// 	return false;
	// }


	private boolean checkAfter(VectorClock prev, String prevOp, ShadowThread currentThread, String curOp, 
			AccessEvent fad, boolean isWrite, ShadowVar p) {

		VectorClock cv = get(currentThread);
		if(prev.anyGt(cv)) { 
			int start=0; 
			while(true) {
				start=prev.nextGt(cv, start);
				if (start==-1) {
					break;
				} 
				Object target = fad.getTarget();
				if (fad.getKind() == Kind.ARRAY) {
					ArrayAccessEvent aae = (ArrayAccessEvent)fad;
					final ArrayAccessInfo arrayAccessInfo = aae.getInfo();
					arrayErrors.error(currentThread, arrayAccessInfo,  
							"Guard State", 	prev, 
							"Array",		Util.objectToIdentityString(target) + "[" + aae.getIndex() +"]", 
							"Locks",		currentThread.getLocksHeld(), 
							"Prev Op",		prevOp+"-by-thread-"+start,  
							"Prev Op CV",	prev, 
							"Cur Op", 		curOp,
							"Cur Op CV", 	cv,
							"Stack",		ShadowThread.stackDumpForErrorMessage(currentThread));
					start++;
					return !arrayErrors.stillLooking(arrayAccessInfo);
				} else {
					FieldInfo fd = ((FieldAccessEvent)fad).getInfo().getField();
					errors.error(currentThread, fd, 
									"Guard State", 	prev, 
									"Class",		target==null?fd.getOwner():target.getClass(), 
									"Field",		Util.objectToIdentityString(target) + "." + fd, 
									"Locks",		currentThread.getLocksHeld(), 
									"Prev Op",		prevOp+"-by-thread-"+start,  
									"Prev Op CV",	prev, 
									"Cur Op", 		curOp,
									"Cur Op CV", 	cv,
									"Stack",		ShadowThread.stackDumpForErrorMessage(currentThread));
					start++;
					return !errors.stillLooking(fd);
				}
			}
			return false;
		} else {
			return false;
		}
	}

	private Integer identitycounter = new Integer(0);

	@Override
	public ShadowVar makeShadowVar(AccessEvent fae) {
		VectorClockPair v = new VectorClockPair();
		synchronized (identitycounter) {
			identitycounter++;
			v.Identity = identitycounter;
		}
		return v;
	}

	@Override
	public void preStart(final StartEvent se) {

		final ShadowThread td = se.getThread();
		final ShadowThread forked = se.getNewThread();

		synchronized (mylog) {
			
			try {
				SigCounter++;
				Integer id = SigCounter;
				mylog.write((td.getTid()+1) + ",SIG," + id + ",nil," + id+"\n");
				mylog.write((forked.getTid()+1) + ",WT," + id + ",nil," + id+"\n");
			} catch (Exception e) {
				System.out.println("bad!");
			}
			//mylog.flush();
		}

		ts_set_cv_hb(forked, new VectorClock(ts_get_cv_hb(td)));

		this.tick(forked);
		this.tick(td);

		super.preStart(se);
	}

	@Override
	public void preNotify(NotifyEvent we) {
		tick(we.getThread());
		// synchronized (mylog) {
		// 	// try {
		// 	// mylog.write((we.getThread().getTid()+1) + ",SIG," + we.getLock().hashCode() + ",nil," + we.getLock().hashCode()+"\n");
		// 	// } catch (Exception e) {
		// 	// System.out.println("bad!");
		// }
			//mylog.flush();
		//}
		synchronized(we.getLock()) {
			get(we.getLock()).max(get(we.getThread()));
		}
		tick(we.getThread());
		super.preNotify(we);
	}

	@Override
	public void preWait(WaitEvent we) {
		tick(we.getThread());
		// synchronized (mylog) {
		// 	// try {
		// 	// mylog.write((we.getThread().getTid()+1) + ",WT," + we.getLock().hashCode() + ",nil," + we.getLock().hashCode()+"\n");
		// 	// } catch (Exception e) {
		// 	// System.out.println("bad!");
		// }
			//mylog.flush();
		//}
		synchronized(we.getLock()) {
			get(we.getLock()).max(get(we.getThread()));
		}
		tick(we.getThread());
		super.preWait(we);
	}

	@Override
	public void postWait(WaitEvent we) { 
		tick(we.getThread());
		synchronized(we.getLock()) {
			get(we.getThread()).max(get(we.getLock()));
		}
		tick(we.getThread());
		super.postWait(we);
	}

	@Override
	public void postJoin(JoinEvent je) { 
		final ShadowThread currentThread = je.getThread();
		final ShadowThread joinedThread = je.getJoiningThread();

		synchronized (mylog) {
			try {
				SigCounter++;
				Integer id = SigCounter;
				mylog.write((currentThread.getTid()+1) + ",WT," + id + ",nil," + id+"\n");
				mylog.write((joinedThread.getTid()+1) + ",SIG," + id + ",nil," + id+"\n");
			} catch (Exception e) {
				System.out.println("bad!");
			}
			//mylog.flush();
		}

		tick(currentThread);
		get(currentThread).max(get(joinedThread));
		tick(currentThread);
		super.postJoin(je);
	}

	private final Decoration<ShadowThread, VectorClock> cvForExit = 
		ShadowThread.makeDecoration("HB:sbarrier", DecorationFactory.Type.MULTIPLE, new NullDefault<ShadowThread, VectorClock>());

	public void postDoBarrier(BarrierEvent<HBBarrierState> be) {
		ShadowThread currentThread = be.getThread();
		VectorClock old = cvForExit.get(currentThread);
		be.getBarrier().reset(old);
		get(currentThread).max(old);
		this.tick(currentThread);
	}

	public void preDoBarrier(BarrierEvent<HBBarrierState> be) {
		ShadowThread td = be.getThread();
		VectorClock entering = be.getBarrier().entering;
		entering.max(get(td));	
		cvForExit.set(td, entering);
	}
	

	@Override
	public ShadowVar cloneState(ShadowVar shadowVar) {
		return null;
	}
	
	
	protected static Decoration<ShadowThread,Vector<VectorClock>> interruptions = 
		ShadowThread.makeDecoration("interruptions", DecorationFactory.Type.MULTIPLE,
				new DefaultValue<ShadowThread, Vector<VectorClock>>() { public Vector<VectorClock> get(ShadowThread ld) { return new Vector<VectorClock>(); }} );

	@Override
	public synchronized void preInterrupt(InterruptEvent me) {
		final ShadowThread td = me.getThread();
		final ShadowThread joining = me.getInterruptedThread();

		VectorClock cv = new VectorClock(get(td));
		interruptions.get(joining).add(cv);
		tick(td);
	}

	@Override
	public synchronized void interrupted(InterruptedEvent e) {
		final ShadowThread current = e.getThread();
		Vector<VectorClock> v = interruptions.get(current);

		for (VectorClock cv : v) {
			get(current).max(cv);
		} 
		
		v.clear();
		super.interrupted(e); 
	}

	@Override
	public void stop(ShadowThread st) {
		synchronized(mylog) {
		try {
			mylog.flush();
			
			// BufferedWriter fileLoc = new BufferedWriter(new FileWriter("files.log"), 32768);
			// Set set = hm.entrySet();
			// Iterator i = set.iterator();
			// while(i.hasNext()) {
			// 	Map.Entry me = (Map.Entry)i.next();
			// 	fileLoc.write(me.getKey() + "," + me.getValue() + "\n");
			// }
			// fileLoc.flush();
			
		} catch (Exception e) {
			System.out.println("bad!");
		}
	}
		//synchronized (maxEpochPerTid) {
		//	maxEpochPerTid.set(st.getTid(), ts_get_E(st));
		//}
		super.stop(st);
		//if (COUNT_OPERATIONS) other.inc(st.getTid());
	}


}
