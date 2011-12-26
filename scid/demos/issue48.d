module scid.demos.issue48;

version( demo ):
import scid.demos.common;

/** Issue 48 - Vector should be a random-access range */
void testIssue48() {
    import std.range;
    static assert(isInputRange!(Vector!double));
    static assert(isForwardRange!(Vector!double));
    static assert(isBidirectionalRange!(Vector!double));
    static assert(isRandomAccessRange!(Vector!double));
}
