package test;

public class Test2 extends Thread {
	static int x;
	static int y;
	static Integer myLock = new Integer(0);

	static void sleepSec(float sec) {
		try{
			Thread.sleep((long)(sec * 1000));
		} catch(InterruptedException e) {
			throw new RuntimeException(e);
		}
	}

	public static class Test3 extends Thread implements Runnable {
		public void run() {
			sleepSec(2);
			synchronized(myLock){
				y = 3;
			}
			x = 3;
		}
	}

	@Override
	public void run() {
		x = 1;
		synchronized(myLock) {
			y = 2;
		}
	}

	public static void main(String args[]) throws Exception {
		final Test2 t1 = new Test2();
		final Test3 t2 = new Test3();

		t1.start();
		t2.start();
		t1.join();
		t2.join();
	}
}
