package main;

import area.Area;
import area.impl.Circle;
import area.impl.Circle2;
import nystrom.NystromMethod;
import nystrom.newnestrom;

public class Main {
    public static void main(String[] args) {
        Area area = new Circle2();
        System.out.println("n\t| u - uExact");
        System.out.println("----------------------------");
        for (int n = 4; n < 200; n = n * 2) {
            NystromMethod nystromMethod = new NystromMethod(n, area, 10);
//            newnestrom nystromMethod = new newnestrom(n, area, 1);
            double uExact = nystromMethod.getExactSolution();
            double u = nystromMethod.getSolution();
            System.out.println(n + "\t| " + Math.abs(u - uExact));
        }

    }
}
