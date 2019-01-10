package test;

public class Test3 extends Thread {
	static int x = 0;
	static int y = 0;
	static int z = 0;

	static Integer myLock = new Integer(0);

	public static class Test4 extends Thread implements Runnable {
		public void run() {
			sleepSec(2);
			synchronized(myLock){
				z = 2;
			}
			x = 2;
		}
	}

	public static class Test5 extends Thread implements Runnable {
		public void run() {
			sleepSec(5);
			z = x; //read x
			y = 3;
		}
	}

	@Override
	public void run() {
		y = 1;
		x = 1;
		synchronized(myLock) {
			z = 1;
		}
	}

	public static void main(String args[]) throws Exception {
		final Test3 t1 = new Test3();
		final Test4 t2 = new Test4();
		final Test5 t3 = new Test5();

		t1.start();
		t2.start();
		t3.start();

		t1.join();
		t2.join();
		t3.join();
	}

	static void sleepSec(float sec) {
		try{
			Thread.sleep((long)(sec * 1000));
		} catch(InterruptedException e) {
			throw new RuntimeException(e);
		}
	}
}


